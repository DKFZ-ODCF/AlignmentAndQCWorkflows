#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

source "$TOOL_WORKFLOW_LIB"

printInfo

set -o pipefail

ON_CONVEY=${useAcceleratedHardware:-false}
ID=${RUN}_${LANE}
SM=sample_${SAMPLE}_${PID}

# RODDY_SCRATCH is used here. It is $PBS_SCRATCH_DIR/$PBS_JOBID for Torque/PBS, and /tmp/roddyScratch/jobid for SGE.
RODDY_BIG_SCRATCH=$(getBigScratchDirectory "${FILENAME_SORTED_BAM}_TEMP")
mkdir -p "$RODDY_BIG_SCRATCH"

# pipes via local scratch dir
FNPIPE1=${RODDY_SCRATCH}/NAMED_PIPE1
FNPIPE2=${RODDY_SCRATCH}/NAMED_PIPE2
NP_READBINS_IN=${RODDY_SCRATCH}/np_readbins_in
NP_COVERAGEQC_IN=${RODDY_SCRATCH}/np_coverageqc_in
NP_COMBINEDANALYSIS_IN=${RODDY_SCRATCH}/np_combinedanalysis_in
NP_FLAGSTATS=${RODDY_SCRATCH}/np_flagstats_in
NP_SAMTOOLS_INDEX_IN=${RODDY_SCRATCH}/np_samtools_index_in

MBUF_SMALL="${MBUFFER_BINARY} -m $MBUFFER_SIZE_SMALL -q -l /dev/null"
MBUF_LARGE="${MBUFFER_BINARY} -m $MBUFFER_SIZE_LARGE -q -l /dev/null"

mkfifo ${NP_READBINS_IN} ${NP_COVERAGEQC_IN} ${NP_COMBINEDANALYSIS_IN} ${NP_FLAGSTATS}

bamname=`basename ${FILENAME_SORTED_BAM}`
INDEX_FILE=${FILENAME_SORTED_BAM}.bai
tempSortedBamFile=${FILENAME_SORTED_BAM}.tmp
tempFileForSort=${RODDY_BIG_SCRATCH}/${bamname}_forsorting
tempBamIndexFile=${FILENAME_SORTED_BAM}.tmp.bai
tempFlagstatsFile=${FILENAME_FLAGSTATS}.tmp

RAW_SEQ=${RAW_SEQ_1}
source ${TOOL_COMMON_ALIGNMENT_SETTINGS_SCRIPT}

NP_BAMSORT=${RODDY_SCRATCH}/NAMED_PIPE_BAMSORT
mkfifo ${NP_BAMSORT}

# Error tracking
FILENAME_BWA_ERRORCODE="$DIR_TEMP/${bamname}_ec_bbam"
FILENAME_BWA_LOG="$DIR_TEMP/${bamname}_errlog_bwamem"
# samtools sort may complain about truncated temp files and for each line outputs
# the error message. This happens when the same files are written at the same time by independent processes,
# see http://sourceforge.net/p/samtools/mailman/samtools-help/thread/BAA90EF6FE3B4D45A7B2F6E0EC5A8366DA3AB5@USTLMLLYC102.rf.lilly.com/
FILENAME_SORT_LOG="$DIR_TEMP/${bamname}_errlog_sort"


bamFileExists=false
# in case the BAM already exists, but QC files are missing, create these only
if [[ -f ${FILENAME_SORTED_BAM} ]] && [[ -s ${FILENAME_SORTED_BAM} ]]
then
    checkBamIsComplete "$FILENAME_SORTED_BAM"
	bamFileExists=true
fi

# test if one of the fastq files is a fake fastq file to simulate paired end sequencing in PIDs with mixed sequencing (single and paired end)
LENGTH_SEQ_1=`${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} 2>/dev/null | head | wc -l`
LENGTH_SEQ_2=`${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} 2>/dev/null | head | wc -l`
[[ ${LENGTH_SEQ_1} == 0 || ${LENGTH_SEQ_2} == 0 ]] && useSingleEndProcessing=true

# make biobambam sort default
useBioBamBamSort="${useBioBamBamSort-true}"
# Do not use adaptor trimming by default.
useAdaptorTrimming="${useAdaptorTrimming-false}"

# Default: Dummy process IDs to simplify downstream logic.
# TODO Remove after completely switching to PID registry system.
true & procUnpack1=$!
true & procUnpack2=$!

fqName="fastq"
mkPipePairSource "$fqName"
if [[ "$bamFileExists" == "false" ]]; then	# we have to make the BAM
    getFastqAsciiStream "$RAW_SEQ_1" > $(getPairedPipeEndPath 1 "$fqName") & procUnpack1=$!
    getFastqAsciiStream "$RAW_SEQ_2" > $(getPairedPipeEndPath 2 "$fqName") & procUnpack2=$!

    if [[ "$qualityScore" == "illumina" ]]; then
        extendPipe $(mkPairedPipeName 1 "$fqName") "qScore" -- toIlluminaScore
        extendPipe $(mkPairedPipeName 2 "$fqName") "qScore" -- toIlluminaScore
    fi

    if [[ "${useAdaptorTrimming:-false}" == "true" ]]; then
        extendPipePair "$fqName" "trimmomatic" -- trimmomatic
    fi

    INPUT_PIPES=""
    if [[ ${LENGTH_SEQ_1} == 0 ]]; then
        cat "$(getPairedPipeEndPath 1 $fqName)" > /dev/null
    else
	    INPUT_PIPES="$(getPairedPipeEndPath 1 $fqName)"
    fi

    if [[ ${LENGTH_SEQ_2} == 0 ]]; then
        cat "$(getPairedPipeEndPath 2 $fqName)" > /dev/null
    else
	    INPUT_PIPES="${INPUT_PIPES} $(getPairedPipeEndPath 2 $fqName)"
    fi

    if [[ -z "$INPUT_PIPES" ]]; then
	    throw 100 "Your input FASTQs are empty: '$RAW_SEQ_1', '$RAW_SEQ_2'"
	fi

fi

# Try to read from pipes BEFORE they are filled.
# In all cases:
# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN} -a ${FILENAME_DIFFCHROM_MATRIX}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX}.tmp  -f ${FILENAME_EXTENDED_FLAGSTATS}.tmp  -m ${FILENAME_ISIZES_STATISTICS}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBA=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN} --outputFile=${FILENAME_GENOME_COVERAGE}.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage=$!

# read bins
(set -o pipefail; ${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE} | $MBUF_SMALL | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}.tmp) & procIDReadbinsCoverage=$!

# use sambamba for flagstats
${SAMBAMBA_FLAGSTATS_BINARY} flagstat -t 1 "$NP_FLAGSTATS" > "$tempFlagstatsFile" & procIDFlagstat=$!

if [[ ${bamFileExists} == true ]]
then
	echo "The BAM file already exists, re-creating other output files."
	# make all the pipes
	(cat ${FILENAME_SORTED_BAM} \
	    | ${MBUF_LARGE} \
	    | tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} ${NP_FLAGSTATS} \
	    | ${SAMBAMBA_BINARY} view /dev/stdin \
	    | ${MBUF_LARGE} > $NP_COMBINEDANALYSIS_IN) \
	    & procIDOutPipe=$!

	wait $procIDOutPipe || throw 13 "Error from sambamba view pipe"

else
	if [[ "$ON_CONVEY" == true ]]
	then	# we have to use sambamba and cannot make an index (because sambamba does not work with a pipe)
		# Here, we always use the local scratch (${RODDY_SCRATCH}) for sorting!
		useBioBamBamSort=false;
		(set -o pipefail; ${BWA_ACCELERATED_BINARY} mem ${BWA_MEM_CONVEY_ADDITIONAL_OPTIONS} \
		    -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" $BWA_MEM_OPTIONS "$INDEX_PREFIX" ${INPUT_PIPES} 2> $FILENAME_BWA_LOG \
		    | $MBUF_LARGE \
		    | tee $NP_COMBINEDANALYSIS_IN \
		    | ${SAMBAMBA_BINARY} view -f bam -S -l 0 -t 8 /dev/stdin \
		    | $MBUF_LARGE \
		    | tee ${NP_FLAGSTATS} \
		    | ${SAMBAMBA_BINARY} sort --tmpdir=${RODDY_SCRATCH} -l 9 -t ${CONVEY_SAMBAMBA_SAMSORT_THREADS} \
		      -m ${CONVEY_SAMBAMBA_SAMSORT_MEMSIZE} /dev/stdin -o /dev/stdout 2> $FILENAME_SORT_LOG \
		    | tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} > $tempSortedBamFile; \
		  echo "$? (Pipe Exit Codes: ${PIPESTATUS[@]})" > "$FILENAME_BWA_ERRORCODE") & procID_MEMSORT=$!

	elif [[ ${useBioBamBamSort} == false ]]
	then	# we use samtools for making the index
	    NP_SORT_ERRLOG="$RODDY_SCRATCH/NP_SORT_ERRLOG"
		mkfifo $NP_SORT_ERRLOG ${NP_SAMTOOLS_INDEX_IN}
		${SAMTOOLS_BINARY} index ${NP_SAMTOOLS_INDEX_IN} ${tempBamIndexFile} & procID_IDX=$!
		(set -o pipefail; ${BWA_BINARY} mem -t ${BWA_MEM_THREADS} \
		    -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" $BWA_MEM_OPTIONS ${INDEX_PREFIX} ${INPUT_PIPES} 2> $FILENAME_BWA_LOG \
		    | $MBUF_LARGE \
		    | tee $NP_COMBINEDANALYSIS_IN \
		    | ${SAMTOOLS_BINARY} view -uSbh - \
		    | $MBUF_LARGE \
		    | ${SAMTOOLS_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} -o - ${tempFileForSort} 2> $NP_SORT_ERRLOG \
		    | tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} ${NP_FLAGSTATS} ${NP_SAMTOOLS_INDEX_IN} > ${tempSortedBamFile}; \
		  echo "$? (Pipe Exit Codes: ${PIPESTATUS[@]})" > "$FILENAME_BWA_ERRORCODE") & procID_MEMSORT=$!

   		# filter samtools error log
		(cat $NP_SORT_ERRLOG | uniq > $FILENAME_SORT_LOG) & procID_logwrite=$!
		wait $procID_logwrite	# do we need a check for it?
		wait $procID_IDX || throw 10 "Error from samtools index"

	else	# biobambam makes the index
		(cat ${NP_BAMSORT} | tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} ${NP_FLAGSTATS} > ${tempSortedBamFile}) & procIDview=$!
		# Output sam to separate named pipe
		# Rewrite to a bamfile
		(set -o pipefail; ${BWA_BINARY} mem -t ${BWA_MEM_THREADS} \
		    -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" $BWA_MEM_OPTIONS ${INDEX_PREFIX} ${INPUT_PIPES} 2> $FILENAME_BWA_LOG \
		    | $MBUF_LARGE \
		    | tee $NP_COMBINEDANALYSIS_IN \
		    | ${SAMTOOLS_BINARY} view -uSbh - \
		    | $MBUF_LARGE \
		    | ${BAMSORT_BINARY} O=${NP_BAMSORT} level=1 inputthreads=2 outputthreads=2 \
		        index=1 indexfilename=${tempBamIndexFile} calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${INDEX_PREFIX} \
		        tmpfile=${tempFileForSort} 2> $FILENAME_SORT_LOG; \
		  echo "$? (Pipe Exit Codes: ${PIPESTATUS[@]})" > "$FILENAME_BWA_ERRORCODE") & procID_MEMSORT=$!


		wait $procIDview || throw 12 "Error from cat from bamsort output for pipes"
	fi

	# Rename BAM file when it has been produced correctly.
    wait $procID_MEMSORT || throw 11 "Error from BWA MEM and sorting pipeline"
	[[ -p $i1 ]] && rm $i1 $i2 $o1 $o2 2> /dev/null
	rm $FNPIPE1
	rm $FNPIPE2
	checkBwaOutput "$tempSortedBamFile" "$FILENAME_BWA_LOG" "$FILENAME_SORT_LOG" "$FILENAME_BWA_ERRORCODE"
	mv "$tempSortedBamFile" "$FILENAME_SORTED_BAM" || throw 36 "Could not move file"
	# The index is only created by samtools or biobambam when producing the BAM, it may be older than the BAM, so update time stamp.
	if [[ -f "$tempBamIndexFile" ]]; then
        mv "$tempBamIndexFile" "$INDEX_FILE" && touch "$INDEX_FILE" || throw 37 "Could not move file"
	fi
	# Clean up adapter trimming pipes if they exist.
	[[ -p $i1 ]] && rm $i1 $i2 $o1 $o2 2> /dev/null
fi

wait $procTrim || throw 38 "Error from trimming"
wait $procIDFlagstat; [[ $? -gt 0 ]] && echo "Error from sambamba flagstats" && exit 14
wait $procIDReadbinsCoverage; [[ $? -gt 0 ]] && echo "Error from genomeCoverage read bins" && exit 15
wait $procIDGenomeCoverage; [[ $? -gt 0 ]] && echo "Error from coverageQCD" && exit 16
wait $procIDCBA; [[ $? -gt 0 ]] && echo "Error from combined QC perl script" && exit 17

# rename QC files
mv ${FILENAME_DIFFCHROM_MATRIX}.tmp ${FILENAME_DIFFCHROM_MATRIX} || throw 28 "Could not move file"
mv ${FILENAME_ISIZES_MATRIX}.tmp ${FILENAME_ISIZES_MATRIX} || throw 29 "Could not move file"
mv ${FILENAME_EXTENDED_FLAGSTATS}.tmp ${FILENAME_EXTENDED_FLAGSTATS} || throw 30 "Could not move file"
mv ${FILENAME_ISIZES_STATISTICS}.tmp ${FILENAME_ISIZES_STATISTICS} || throw 31 "Could not move file"
mv ${FILENAME_DIFFCHROM_STATISTICS}.tmp ${FILENAME_DIFFCHROM_STATISTICS} || throw 32 "Could not move file"
mv ${FILENAME_READBINS_COVERAGE}.tmp ${FILENAME_READBINS_COVERAGE} || throw 34 "Could not move file"
mv ${FILENAME_GENOME_COVERAGE}.tmp ${FILENAME_GENOME_COVERAGE} || throw 35 "Could not move file"
mv ${tempFlagstatsFile} ${FILENAME_FLAGSTATS} || throw 33 "Could not move file"

runFingerprinting "${FILENAME_SORTED_BAM}" "$FILENAME_FINGERPRINTS"
removeRoddyBigScratch

# QC summary
# remove old warnings file if it exists (due to errors in run such as wrong chromsizes file)
[[ -f ${FILENAME_QCSUMMARY}_WARNINGS.txt ]] && rm ${FILENAME_QCSUMMARY}_WARNINGS.txt

${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY \
    -p "$PID" \
    -s "$SAMPLE" \
    -r "$RUN" \
    -l "$LANE" \
    -w "${FILENAME_QCSUMMARY}_WARNINGS.txt" \
    -f "$FILENAME_FLAGSTATS" \
    -d "$FILENAME_DIFFCHROM_STATISTICS" \
    -i "$FILENAME_ISIZES_STATISTICS" \
    -c "$FILENAME_GENOME_COVERAGE" \
    > "${FILENAME_QCSUMMARY}_temp" \
    || throw 14 "Error from writeQCsummary.pl"
mv "${FILENAME_QCSUMMARY}_temp" "$FILENAME_QCSUMMARY"

groupLongAndShortChromosomeNames "$FILENAME_GENOME_COVERAGE" \
    > "$FILENAME_GROUPED_GENOME_COVERAGE.tmp"  \
    || throw 43 "Error grouping reads by having (=long) or not having (=short) prefix/suffix"
mv "$FILENAME_GROUPED_GENOME_COVERAGE.tmp" "$FILENAME_GROUPED_GENOME_COVERAGE" \
    || throw 27 "Could not move file"

# Produce qualitycontrol.json for OTP.
${PERL_BINARY} ${TOOL_QC_JSON} \
    "$FILENAME_GENOME_COVERAGE" \
    "$FILENAME_GROUPED_GENOME_COVERAGE" \
    "$FILENAME_ISIZES_STATISTICS" \
    "$FILENAME_FLAGSTATS" \
    "$FILENAME_DIFFCHROM_STATISTICS" \
    > "$FILENAME_QCJSON.tmp" \
    || throw 25 "Error when compiling qualitycontrol.json '${FILENAME_QCJSON}', stopping here"

mv ${FILENAME_QCJSON}.tmp ${FILENAME_QCJSON} || throw 27 "Could not move file"

# plots are only made for paired end and not on convey
if [[ ${useSingleEndProcessing-false} == true ]] || [[ "$ON_CONVEY" == "true" ]]; then
    exit 0
else

    $RSCRIPT_BINARY $TOOL_INSERT_SIZE_PLOT_SCRIPT \
        "$FILENAME_ISIZES_MATRIX" \
        "$FILENAME_ISIZES_STATISTICS" \
        "${FILENAME_ISIZES_PLOT}_temp" \
        "PE insertsize of $bamname" \
        || throw 22 "Error from insert sizes plotter"

    mv "${FILENAME_ISIZES_PLOT}_temp" "$FILENAME_ISIZES_PLOT" || throw 27 "Could not move file"

    $RSCRIPT_BINARY $TOOL_PLOT_DIFFCHROM \
        -i "$FILENAME_DIFFCHROM_MATRIX" \
        -s "$FILENAME_DIFFCHROM_STATISTICS" \
        -o "${FILENAME_DIFFCHROM_PLOT}_temp" \
        || throw 23 "Error from chrom_diff.r"

    mv  "${FILENAME_DIFFCHROM_PLOT}_temp" "$FILENAME_DIFFCHROM_PLOT" || throw 27 "Could not move file"

fi