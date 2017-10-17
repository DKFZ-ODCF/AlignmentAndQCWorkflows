#!/bin/bash

# Slim version of sampe sort

set -o pipefail

ID=${RUN}_${LANE}
SM=sample_${SAMPLE}_${PID}

# use scratch dir for temp files: samtools sort uses the current working directory for them unless provided with a different prefix
# so far, it's not the local scratch (${RODDY_SCRATCH}) but inside roddyExecutionStore of the PID
# here, using the job ID is a bad idea because if the job fails, the temp files remain
#WORKDIR=${DIR_TEMP}/${RODDY_JOBID}
WORKDIR=${FILENAME_SORTED_BAM}_TEMP_SORTING
mkdir -p ${WORKDIR}

# pipes via local scratch dir
FNPIPE1=$RODDY_SCRATCH/NAMED_PIPE1
FNPIPE2=$RODDY_SCRATCH/NAMED_PIPE2
NP_READBINS_IN=${RODDY_SCRATCH}/np_readbins_in
NP_COVERAGEQC_IN=${RODDY_SCRATCH}/np_coverageqc_in
NP_COMBINEDANALYSIS_IN=${RODDY_SCRATCH}/np_combinedanalysis_in
NP_FLAGSTATS_IN=${RODDY_SCRATCH}/np_flagstats_in
NP_INDEX_IN=${RODDY_SCRATCH}/np_index_in

MBUF_100M="${MBUFFER_BINARY} -m 100m -q -l /dev/null"
MBUF_1G="${MBUFFER_BINARY} -m 1G -q -l /dev/null"
MBUF_2G="${MBUFFER_BINARY} -m 2G -q -l /dev/null"
MBUF_4G="${MBUFFER_BINARY} -m 4G -q -l /dev/null"

mkfifo ${NP_READBINS_IN} ${NP_COVERAGEQC_IN} ${NP_COMBINEDANALYSIS_IN} ${NP_FLAGSTATS_IN} ${NP_INDEX_IN}

bamname=`basename ${FILENAME_SORTED_BAM}`
FILENAME_BWA_LOG=${DIR_TEMP}/${bamname}_errlog_sampe
# logfile for samtools sort, may complain about truncated temp files and for each line outputs
# the error message. This happens when the same files are written at the same time,
# see http://sourceforge.net/p/samtools/mailman/samtools-help/thread/BAA90EF6FE3B4D45A7B2F6E0EC5A8366DA3AB5@USTLMLLYC102.rf.lilly.com/
# To redirect STDERR and keep the log file small,
# use named pipe as outfile and run uniq over it
NP_SORT_ERRLOG=$RODDY_SCRATCH/NP_SORT_ERRLOG
mkfifo $NP_SORT_ERRLOG
FILENAME_SORT_LOG=${DIR_TEMP}/${bamname}_errlog_sort

#Samtools always attaches .bam, so we use a different output filename for samtools without .bam
bamIndexFile=${FILENAME_SORTED_BAM}.bai
tempBamIndexFile=${FILENAME_SORTED_BAM}.tmp.bai
tempSortedBamFile=${FILENAME_SORTED_BAM}_tmp
tempFileForSort=${WORKDIR}/sorted.temporary
TMP_FILE=${tempSortedBamFile}       # This variable is used for the error checking code.

RAW_SEQ=${RAW_SEQ_1}
source ${TOOL_COMMON_ALIGNMENT_SETTINGS_SCRIPT}

bamFileExists=false
# in case the BAM already exists, but QC files are missing, create these only
if [[ -f ${FILENAME_SORTED_BAM} ]] && [[ -s ${FILENAME_SORTED_BAM} ]]
then
	bamFileExists=true
fi

# Try to read from pipes BEFORE they are filled.
# in all cases:
# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN} -a ${FILENAME_DIFFCHROM_MATRIX}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX}.tmp  -f ${FILENAME_EXTENDED_FLAGSTATS}.tmp  -m ${FILENAME_ISIZES_STATISTICS}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBA=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN} --outputFile=${FILENAME_GENOME_COVERAGE}.tmp --processors=4 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage=$!

# read bins
(set -o pipefail; ${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN} --outputFile=/dev/stdout --processors=2 --mode=countReads --windowSize=${WINDOW_SIZE} | $MBUF_100M | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}.tmp) & procIDReadbinsCoverage=$!

# use sambamba for flagstats
(set -o pipefail; cat ${NP_FLAGSTATS_IN} | ${SAMBAMBA_BINARY} flagstat /dev/stdin > ${FILENAME_FLAGSTATS}.tmp) & procIDFlagstat=$!

# samtools for index
${SAMTOOLS_BINARY} index $NP_INDEX_IN $tempBamIndexFile & procIDindexing=$!

useSingleEndProcessing=${useSingleEndProcessing-false}
if [[ ${bamFileExists} == true ]]
then
	echo "the BAM file already exists, re-creating other output files."
	# make all the pipes
	(cat ${FILENAME_SORTED_BAM} | ${MBUF_2G} | tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} ${NP_FLAGSTATS_IN} | ${MBUF_2G} | ${SAMBAMBA_BINARY} view /dev/stdin | ${MBUF_2G} > $NP_COMBINEDANALYSIS_IN) & procIDOutPipe=$!

else	# we have to make the BAM
	mkfifo $FNPIPE1 $FNPIPE2
	if [[ ${useSingleEndProcessing} == "true" ]]
	then
		# Single end read is active
		# are ${TRIM_STEP} and ${REVERSE_STEP} set correctly?
		nice ${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_1 ${TRIM_STEP} ${REVERSE_STEP} | $MBUF_1G > $FNPIPE1 & procID_singleendin=$!

		(set -o pipefail; ${BWA_BINARY} samse -T -t 8 -r "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" \
		${INDEX_PREFIX} ${FILENAME_SEQ_1} $FNPIPE1 2> ${FILENAME_BWA_LOG} | $MBUF_2G | \
		tee $NP_COMBINEDANALYSIS_IN | ${SAMTOOLS_BINARY} view -uSbh - | $MBUF_2G | \
		tee $NP_FLAGSTATS_IN | ${SAMTOOLS_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} -o - ${tempFileForSort} 2>$NP_SORT_ERRLOG | \
		tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} $NP_INDEX_IN > ${tempSortedBamFile}; echo $? > ${DIR_TEMP}/${bamname}_ec) & procID_SAMSESORT=$!
		# filter samtools error log
		(cat $NP_SORT_ERRLOG | uniq > $FILENAME_SORT_LOG) & procID_logwrite=$!
		wait $procID_logwrite	# do we need a check for it?
		wait $procID_singleendin; [[ $? -gt 0 ]] && echo "Error from single end preprocessing" && exit 10
		wait $procID_SAMSESORT; [[ ! `cat ${DIR_TEMP}/${bamname}_ec` -eq "0" ]] && echo "bwa samse - samtools sort pipe returned a non-zero exit code and the job will die now." && exit 100
	else
		if [[ $useAdaptorTrimming == false ]]
		then
			##Reset error and tmp variables. these are modified by bwaCommonAlignment...
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_1 ${TRIM_STEP} ${REVERSE_STEP} | $MBUF_1G > ${FNPIPE1}" &

			# Repeat some steps for the second named pipe
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_2 ${TRIM_STEP} ${REVERSE_STEP} | $MBUF_1G > ${FNPIPE2}" &

		else
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | $MBUF_4G > $i1" &
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | $MBUF_4G > $i2" &

			#We need a second output pipe for sampe (default is to go to /dev/null from aln)
			o2=${DIR_SCRATCH}/at_o2
			mkfifo $o2

			if [[ "$ADAPTOR_TRIMMING_TOOL" == *.jar ]]
			then
			eval "java7 -jar  ${TOOL_ADAPTOR_TRIMMING} $ADAPTOR_TRIMMING_OPTIONS_0 $i1 $i2 $o1 $u1 $o2 $u2 $ADAPTOR_TRIMMING_OPTIONS_1" &
			fi

			cat $o1 ${TRIM_STEP} ${REVERSE_STEP} | $MBUF_4G > $FNPIPE1 &
			cat $o2 ${TRIM_STEP} ${REVERSE_STEP} | $MBUF_4G > $FNPIPE2 &

		fi
		(set -o pipefail; ${BWA_BINARY} sampe -P -T -t 8 ${BWA_SAMPESORT_OPTIONS} -r "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" \
			${INDEX_PREFIX} ${FILENAME_SEQ_1} ${FILENAME_SEQ_2} $FNPIPE1 $FNPIPE2 2> ${FILENAME_BWA_LOG} | \
			$MBUF_2G | tee $NP_COMBINEDANALYSIS_IN | ${SAMTOOLS_BINARY} view -uSbh - | $MBUF_2G | \
			tee $NP_FLAGSTATS_IN | ${SAMTOOLS_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} -o - ${tempFileForSort} 2>$NP_SORT_ERRLOG | \
			tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} $NP_INDEX_IN > ${tempSortedBamFile}; echo $? > ${DIR_TEMP}/${bamname}_ec) & procID_SAMPESORT=$!
		# filter samtools error log
		(cat $NP_SORT_ERRLOG | uniq > $FILENAME_SORT_LOG) & procID_logwrite=$!
		wait $procID_logwrite	# do we need a check for it?
		wait $procID_SAMPESORT; [[ ! `cat ${DIR_TEMP}/${bamname}_ec` -eq "0" ]] && echo "bwa sampe - samtools sort pipe returned a non-zero exit code and the job will die now." && exit 100
	fi
fi

#Wait for unfinished samtools jobs.
if [[ ${bamFileExists} == false ]]	# BAM was created
then
	rm $FNPIPE1
	rm $FNPIPE2
	errorString="There was a non-zero exit code in the bwa sampe - samtools sort pipeline; exiting..." 
	source ${TOOL_BWA_ERROR_CHECKING_SCRIPT}
	mv ${tempSortedBamFile} ${FILENAME_SORTED_BAM}
fi

#rm $FILE_BENCHMARK_STAYALIVE
#wait $utilizationLogProcess
rm -rf ${DIR_TEMP}/${bamname}

wait $procIDFlagstat; [[ $? -gt 0 ]] && echo "Error from sambamba flagstats" && exit 11
wait $procIDReadbinsCoverage; [[ $? -gt 0 ]] && echo "Error from genomceCoverage read bins" && exit 12
wait $procIDGenomeCoverage; [[ $? -gt 0 ]] && echo "Error from coverageQCD" && exit 13
wait $procIDCBA; [[ $? -gt 0 ]] && echo "Error from combined QC perl script" && exit 14
wait $procIDindexing; [[ $? -gt 0 ]] && echo "Error from samtools index" && exit 16

# update index time stamp
mv ${tempBamIndexFile} ${bamIndexFile} && touch ${bamIndexFile}

# rename QC files
mv ${FILENAME_DIFFCHROM_MATRIX}.tmp ${FILENAME_DIFFCHROM_MATRIX}
mv ${FILENAME_ISIZES_MATRIX}.tmp ${FILENAME_ISIZES_MATRIX}
mv ${FILENAME_EXTENDED_FLAGSTATS}.tmp ${FILENAME_EXTENDED_FLAGSTATS}
mv ${FILENAME_ISIZES_STATISTICS}.tmp ${FILENAME_ISIZES_STATISTICS}
mv ${FILENAME_DIFFCHROM_STATISTICS}.tmp ${FILENAME_DIFFCHROM_STATISTICS}
mv ${FILENAME_READBINS_COVERAGE}.tmp ${FILENAME_READBINS_COVERAGE}
mv ${FILENAME_GENOME_COVERAGE}.tmp ${FILENAME_GENOME_COVERAGE}
mv ${FILENAME_FLAGSTATS}.tmp ${FILENAME_FLAGSTATS}

# QC summary
# remove old warnings file if it exists (due to errors in run such as wrong chromsizes file)
[[ -f ${FILENAME_QCSUMMARY}_WARNINGS.txt ]] && rm ${FILENAME_QCSUMMARY}_WARNINGS.txt
(${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY -p $PID -s $SAMPLE -r $RUN -l $LANE -w ${FILENAME_QCSUMMARY}_WARNINGS.txt -f $FILENAME_FLAGSTATS -d $FILENAME_DIFFCHROM_STATISTICS -i $FILENAME_ISIZES_STATISTICS -c $FILENAME_GENOME_COVERAGE > ${FILENAME_QCSUMMARY}_temp && mv ${FILENAME_QCSUMMARY}_temp $FILENAME_QCSUMMARY) || ( echo "Error from writeQCsummary.pl" && exit 17)

# plots are only made for paired end
[[ ${useSingleEndProcessing-false} == true ]] && exit 0

${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAME_ISIZES_MATRIX} ${FILENAME_ISIZES_STATISTICS} ${FILENAME_ISIZES_PLOT}_temp "PE insertsize of ${bamname}" && mv ${FILENAME_ISIZES_PLOT}_temp ${FILENAME_ISIZES_PLOT} || ( echo "Error from insert sizes plotter" && exit 18)

${RSCRIPT_BINARY} ${TOOL_PLOT_DIFFCHROM} -i "$FILENAME_DIFFCHROM_MATRIX" -s "$FILENAME_DIFFCHROM_STATISTICS" -o "${FILENAME_DIFFCHROM_PLOT}_temp" && mv  ${FILENAME_DIFFCHROM_PLOT}_temp ${FILENAME_DIFFCHROM_PLOT} || ( echo "Error from chrom_diff.r" && exit 19)
