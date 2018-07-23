#!/bin/bash

source "$TOOL_BASH_LIB"
printInfo

set -o pipefail

today=`date +'%Y-%m-%d_%Hh%M'`

localScratchDirectory=${RODDY_SCRATCH}	# is for PBS $PBS_SCRATCHDIR/$PBS_JOBID, for SGE /tmp/roddyScratch/jobid
tempDirectory=${FILENAME}_MOMDUP

tempBamFile=${tempDirectory}/$(basename "$FILENAME")_TMP
tempFlagstatsFile=${FILENAME_FLAGSTATS}.tmp
tempIndexFile=${tempBamFile}.bai.tmp
tempMd5File=${tempBamFile}.md5.tmp

CHECKSUM_BINARY=${CHECKSUM_BINARY-md5sum}

NP_PIC_OUT=${localScratchDirectory}/np_picard_out.sam
NP_SAM_IN=${localScratchDirectory}/np_compression_in
NP_INDEX_IN=${localScratchDirectory}/np_index_in
NP_FLAGSTATS_IN=${localScratchDirectory}/np_flagstats_in
NP_READBINS_IN=${localScratchDirectory}/np_readbins_in
NP_COVERAGEQC_IN=${localScratchDirectory}/np_coverageqc_in
NP_COMBINEDANALYSIS_IN=${localScratchDirectory}/np_combinedanalysis_in
NP_MD5_IN=${localScratchDirectory}/np_md5_in

returnCodeMarkDuplicatesFile=${tempDirectory}/rcMarkDup.txt
bamname=`basename ${FILENAME}`

declare -a INPUT_FILES="$INPUT_FILES"
# Handle existing BAM provided by 'bam' parameter or present as target FILENAME.
if [[ -v bam && ! -z "$bam" ]]; then
    checkBamIsComplete "$bam"
    EXISTING_BAM="$bam"
fi
# Handle existing merged BAM at target location.
if [[ -f ${FILENAME} && -s ${FILENAME} ]]; then
    checkBamIsComplete "$FILENAME"
    if [[ -v EXISTING_BAM ]]; then
        echo "Input BAM provided via 'bam' option has precedence over existing merged BAM. Rescuing merged BAM."
        mv "$FILENAME" "${FILENAME}_before_${today}" || throw 50 "Could not move $FILENAME"
    else
        # If the merged file already exists, only merge new lanes to it.
        EXISTING_BAM="$FILENAME"
    fi
fi

if [[ -v EXISTING_BAM ]]; then
        singlebams=`${SAMTOOLS_BINARY} view -H ${EXISTING_BAM} | grep "^@RG"`
        [[ -z "$singlebams" ]] && echo "could not detect single lane BAM files in ${EXISTING_BAM}, stopping here" && exit 23

        notyetmerged=`perl ${TOOL_CHECK_ALREADY_MERGED_LANES} $(stringJoin ":" ${INPUT_FILES[@]}) "$singlebams" ${pairedBamSuffix} $SAMPLE`
        [[ "$?" != 0 ]] && echo "something went wrong with the detection of merged files in ${EXISTING_BAM}, stopping here" && exit 24

        # the Perl script returns BAM names separated by :, ending with :
        if [ -z $notyetmerged ]
        then
            echo "All listed BAM files are already in ${EXISTING_BAM}, re-creating other output files."
            bamFileExists=true
        else
            echo "New lane(s) need to be merged to the BAM."
		    mkdir -p $tempDirectory
            # input files is now the merged file and the new file(s)
            declare -a INPUT_FILES=("$EXISTING_BAM" $(echo $notyetmerged | sed -re 's/:/ /g'))
            # keep the old metrics file for comparison. For an externally provided BAM (via cvalue 'bam')
            # that file may not exist in the target directory.
            if [[ -f $FILENAME_METRICS ]]; then
                mv ${FILENAME_METRICS} ${FILENAME_METRICS}_before_${today}.txt || throw 37 "Could not move file"
            fi
        fi
else	# BAM needs to be created de novo
	mkdir -p $tempDirectory
fi


mkfifo ${NP_PIC_OUT} ${NP_SAM_IN} ${NP_INDEX_IN} ${NP_FLAGSTATS_IN} ${NP_READBINS_IN} ${NP_COVERAGEQC_IN} ${NP_COMBINEDANALYSIS_IN} ${NP_MD5_IN}

mergeINPUT=''

MBUF_100M="${MBUFFER_BINARY} -m 100m -q -l /dev/null"
MBUF_2G="${MBUFFER_BINARY} -m 2g -q -l /dev/null"
MBUF_4G="${MBUFFER_BINARY} -m 4g -q -l /dev/null"

# default would be to use biobambam
useBioBamBamMarkDuplicates=${useBioBamBamMarkDuplicates-true}

bamFileExists=${bamFileExists-false}

if [[ ${useBioBamBamMarkDuplicates} == false ]]; then
    PICARD_BINARY=${PICARD_BINARY-picard-1.125.sh}

    # The filehandles are used for some picard optimization.
    FILEHANDLES=$((`ulimit -n` - 16))

    # Complicate error catching because:
    # In some cases picard exits before cat starts to work on the named pipes. If so, the cat processes block and wait for input data forever.
    # TODO Think hard if it helps to just move picard two steps back and to start the cat processes earlier.

    # index piped BAM (always done)
    (${SAMTOOLS_BINARY} index ${NP_INDEX_IN} $tempIndexFile) & procIDSamtoolsIndex=$!

    # make picard write SAM output that is later compressed more efficiently with samtools
    if [[ ${bamFileExists} == false ]]; then
        (JAVA_OPTIONS="-Xms64G -Xmx64G" $PICARD_BINARY MarkDuplicates $(toMinusIEqualsList ${INPUT_FILES[@]}) OUTPUT=${NP_PIC_OUT} METRICS_FILE=${FILENAME_METRICS}.tmp MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${FILEHANDLES} TMP_DIR=${tempDirectory} COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT ${mergeAndRemoveDuplicates_argumentList} ASSUME_SORTED=TRUE CREATE_INDEX=FALSE MAX_RECORDS_IN_RAM=12500000; echo $? > ${returnCodeMarkDuplicatesFile}) & procIDPicard=$!

        # provide named pipes of SAM type
        (cat ${NP_PIC_OUT} | ${MBUF_2G} | tee ${NP_SAM_IN} ${NP_COMBINEDANALYSIS_IN} > /dev/null) & procIDPicardOutPipe=$!

        # calculate MD5 sum
        cat ${NP_MD5_IN} | ${CHECKSUM_BINARY} | cut -d ' ' -f 1 > ${tempMd5File} & procIDMd5=$!

        # convert SAM to BAM
        (${SAMTOOLS_BINARY} view -S -@ 8 -b ${NP_SAM_IN} \
            | ${MBUF_2G} \
            | tee ${NP_INDEX_IN} ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} ${NP_MD5_IN} \
            > ${tempBamFile}) & procIDSamtoolsView=$!
    else
        # To prevent abundancy of ifs, reuse the process id another time.
        (cat ${FILENAME} | ${MBUF_2G} | tee ${NP_INDEX_IN} ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} | ${SAMTOOLS_BINARY} view - > ${NP_COMBINEDANALYSIS_IN}) & procIDSamtoolsView=$!; procIDPicardOutPipe=$procIDSamtoolsView
        procIDMd5=$procIDSamtoolsView
    fi

else
    # biobambam outputs BAM, this has to be converted to SAM for Perl script
    #        rewritebam=1 \
    #        rewritebamlevel=1 \
    # using $workDirectory was a bad idea for biobambam: each time it crashes, there are large files left over
    # and they will never be deleted because the directory is different for another job ID - because do not use the scratch
    # so use $tempDirectory instead
    if [[ ${bamFileExists} == false ]]; then
        (${BAMMARKDUPLICATES_BINARY} \
            M=${FILENAME_METRICS}.tmp \
            tmpfile=${tempDirectory}/biobambammerge.tmp \
            markthreads=8 \
            level=9 \
            index=1 \
            indexfilename=$tempIndexFile \
            $(toMinusIEqualsList ${INPUT_FILES[@]}) \
            O=${NP_PIC_OUT}; echo $? > ${returnCodeMarkDuplicatesFile}) & procIDBBB=$!

        # TODO: The NP readers processes need to be started *before* the pipe-writing process. ^^ Move up!

        # create BAM pipes for flagstat and the two D tools, write to BAM
        # REUSE! NP_SAM_IN for bam to sam conversion
        # was: MBUF_4G, why?
        cat ${NP_PIC_OUT} | ${MBUF_2G} | tee ${NP_SAM_IN} ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} ${NP_MD5_IN} \
            > ${tempBamFile} & procIDBBBOutPipe=$!

        # calculate MD5 sum
        cat ${NP_MD5_IN} | ${CHECKSUM_BINARY} | cut -d ' ' -f 1 > ${tempMd5File} & procIDMd5=$!

        # make a SAM pipe for the Perl tool
        ${SAMTOOLS_BINARY} view ${NP_SAM_IN} | ${MBUF_2G} > ${NP_COMBINEDANALYSIS_IN} & procIDSAMpipe=$!
    else
        (cat ${FILENAME} | ${MBUF_2G} | tee ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} | ${SAMTOOLS_BINARY} view - > ${NP_COMBINEDANALYSIS_IN}) & procIDSAMpipe=$!; procIDBBB=$procIDSAMpipe
        procIDMd5=$procIDSAMpipe
    fi
fi

# in all cases:
# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN} -a ${FILENAME_DIFFCHROM_MATRIX}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX}.tmp  -f ${FILENAME_EXTENDED_FLAGSTATS}.tmp  -m ${FILENAME_ISIZES_STATISTICS}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBA=$!

# make flagstats of piped BAM
(${SAMTOOLS_BINARY} flagstat ${NP_FLAGSTATS_IN} > ${tempFlagstatsFile}) & procIDSamtoolsFlagstat=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN} --outputFile=${FILENAME_GENOME_COVERAGE}.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage=$!

# readbins coverage (1 or 10 kb bins)
# TODO: compress with bgzip, index with tabix, make R script use zipped data
#BGZIP_BINARY=bgzip
#TABIX_BINARY=tabix

#--processors=1
# this part often fails with broken pipe, ?? where this comes from. The mbuffer did not help, maybe --processors=4 does?
(${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE} | $MBUF_100M | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}.tmp) & procIDReadbinsCoverage=$!

# Some waits for parallel processes. This also depends on the used merge binary.

# Picard
if [[ ${useBioBamBamMarkDuplicates} == false ]]; then
    if [[ ${bamFileExists} == false ]]; then
	    wait $procIDPicard
	    [[ ! `cat ${returnCodeMarkDuplicatesFile}` -eq "0" ]] && echo "Picard returned an exit code and the job will die now." && exit 100
	    checkBamIsComplete "$tempBamFile"
        mv ${tempBamFile} ${FILENAME} || throw 38 "Could not move file"
        mv ${FILENAME_METRICS}.tmp ${FILENAME_METRICS} || throw 39 "Could not move file"
        mv ${tempMd5File} ${FILENAME}.md5 || throw 36 "Could not move file"
    fi
	# index is always made again, needs to be updated to not be older than BAM
	# this file does not exist!!!!!!!
	wait $procIDPicardOutPipe ; [[ $? -gt 0 ]] && echo "Error from picard SAM pipe" && exit 5
	wait $procIDSamtoolsIndex ; [[ $? -gt 0 ]] && echo "Error from samtools index pipe" && exit 6
	wait $procIDSamtoolsView ; [[ $? -gt 0 ]] && echo "Error from samtools view pipe" && exit 7
	mv $tempIndexFile ${FILENAME}.bai && touch ${FILENAME}.bai
else
  # Biobambam
    wait $procIDBBB; [[ $? -gt 0 ]] && echo "Error from Biobambam markduplicates" && exit 10
    if [[ ${bamFileExists} == false ]]; then
    	[[ ! `cat ${returnCodeMarkDuplicatesFile}` -eq "0" ]] && echo "Biobambam returned an exit code and the job will die now." && exit 100
    	# always rename BAM, even if other jobs might have failed
    	checkBamIsComplete "$tempBamFile"
    	mv ${tempBamFile} ${FILENAME} || throw 40 "Could not move file"
    	mv $tempIndexFile ${FILENAME}.bai && touch ${FILENAME}.bai || throw 41 "Could not move file"   # Update timestamp because by piping the index might be older than the BAM
        mv ${FILENAME_METRICS}.tmp ${FILENAME_METRICS} || throw 42 "Could not move file"
        mv ${tempMd5File} ${FILENAME}.md5 || throw 36 "Could not move file"
        wait $procIDBBBOutPipe; [[ $? -gt 0 ]] && echo "Error from biobambam BAM pipe" && exit 8
    fi
    wait $procIDSAMpipe; [[ $? -gt 0 ]] && echo "Error from samtools SAM pipe" && exit 9
fi

wait $procIDReadbinsCoverage; [[ $? -gt 0 ]] && echo "Error from genomeCoverage read bins" && exit 10
wait $procIDGenomeCoverage; [[ $? -gt 0 ]] && echo "Error from coverageQCD" && exit 11
wait $procIDSamtoolsFlagstat; [[ $? -gt 0 ]] && echo "Error from samtools flagstats" && exit 12
wait $procIDCBA; [[ $? -gt 0 ]] && echo "Error from combined QC perl script" && exit 13
wait $procIDMd5; [[ $? -gt 0 ]] && echo "Error from ${CHECKSUM_BINARY}" && exit 15

mv ${FILENAME_DIFFCHROM_MATRIX}.tmp ${FILENAME_DIFFCHROM_MATRIX} || throw 28 "Could not move file"
mv ${FILENAME_ISIZES_MATRIX}.tmp ${FILENAME_ISIZES_MATRIX} || throw 29 "Could not move file"
mv ${FILENAME_EXTENDED_FLAGSTATS}.tmp ${FILENAME_EXTENDED_FLAGSTATS} || throw 30 "Could not move file"
mv ${FILENAME_ISIZES_STATISTICS}.tmp ${FILENAME_ISIZES_STATISTICS} || throw 31 "Could not move file"
mv ${FILENAME_DIFFCHROM_STATISTICS}.tmp ${FILENAME_DIFFCHROM_STATISTICS} || throw 32 "Could not move file"
mv ${tempFlagstatsFile} ${FILENAME_FLAGSTATS} || throw 33 "Could not move file"
mv ${FILENAME_READBINS_COVERAGE}.tmp ${FILENAME_READBINS_COVERAGE} || throw 34 "Could not move file"
mv ${FILENAME_GENOME_COVERAGE}.tmp ${FILENAME_GENOME_COVERAGE} || throw 35 "Could not move file"

# QC summary
runExomeAnalysis=${runExomeAnalysis-false}
if [ "$runExomeAnalysis" = "true" ]
then
    analysis_type=exome
else
    analysis_type=genome
fi

# if the warnings file had been created before, remove it:
# it may be from one lane WGS with < 30x (which raises a warning)

[[ -f ${FILENAME_QCSUMMARY}_WARNINGS.txt ]] && rm ${FILENAME_QCSUMMARY}_WARNINGS.txt

if [[ "$bamFileExists" == true || (-v mergeOnly && "$mergeOnly" == true) ]]; then
    METRICS_OPTION="";
else
    METRICS_OPTION="-m ${FILENAME_METRICS}"
fi
${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY -p $PID -s $SAMPLE -r all_merged -l $analysis_type -w ${FILENAME_QCSUMMARY}_WARNINGS.txt -f $FILENAME_FLAGSTATS -d $FILENAME_DIFFCHROM_STATISTICS -i $FILENAME_ISIZES_STATISTICS -c $FILENAME_GENOME_COVERAGE ${METRICS_OPTION} > ${FILENAME_QCSUMMARY}_temp && mv ${FILENAME_QCSUMMARY}_temp $FILENAME_QCSUMMARY || { echo "Error from writeQCsummary.pl" && exit 14; }

[[ -d $tempDirectory ]] && rm -rf $tempDirectory

groupLongAndShortChromosomeNames "$FILENAME_GENOME_COVERAGE" \
    > "$FILENAME_GROUPED_GENOME_COVERAGE.tmp"  \
    || throw 43 "Error grouping reads by having (=long) or not having (=short) prefix/suffix"
mv "$FILENAME_GROUPED_GENOME_COVERAGE.tmp" "$FILENAME_GROUPED_GENOME_COVERAGE" || throw 27 "Could not move file"

# Produced qualitycontrol.json for OTP. Remove the first branch as soon as the dip-statistics is implemented.
${PERL_BINARY} ${TOOL_QC_JSON} \
    ${FILENAME_GENOME_COVERAGE} \
    ${FILENAME_GROUPED_GENOME_COVERAGE} \
    ${FILENAME_ISIZES_STATISTICS} \
    ${FILENAME_FLAGSTATS} \
    ${FILENAME_DIFFCHROM_STATISTICS} \
    > ${FILENAME_QCJSON}.tmp \
    || throw 25 "Error when compiling qualitycontrol.json for ${FILENAME}"
mv ${FILENAME_QCJSON}.tmp ${FILENAME_QCJSON} || throw 27 "Could not move file"

# if the BAM only contains single end reads, there can be no pairs to have insert sizes or ends mapping to different chromosomes
# and plotting would fail
grep -w "NA" $FILENAME_DIFFCHROM_STATISTICS && useSingleEndProcessing=true

[[ ${useSingleEndProcessing-false} == true ]] && exit 0

${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAME_ISIZES_MATRIX} ${FILENAME_ISIZES_STATISTICS} ${FILENAME_ISIZES_PLOT}_temp "PE insertsize of ${bamname} (rmdup)" && mv  ${FILENAME_ISIZES_PLOT}_temp ${FILENAME_ISIZES_PLOT} || throw 22 "Error from insert sizes plotter"
${RSCRIPT_BINARY} ${TOOL_PLOT_DIFFCHROM} -i "$FILENAME_DIFFCHROM_MATRIX" -s "$FILENAME_DIFFCHROM_STATISTICS" -o "${FILENAME_DIFFCHROM_PLOT}_temp" && mv  ${FILENAME_DIFFCHROM_PLOT}_temp ${FILENAME_DIFFCHROM_PLOT} || throw 23 "Error from chrom_diff.r"

