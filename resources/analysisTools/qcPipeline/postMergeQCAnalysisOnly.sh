#!/bin/bash

set -o pipefail
set -u
#~ today=`date +'%Y-%m-%d_%Hh%M'`
localScratchDirectory=${RODDY_SCRATCH}	# is for PBS $PBS_SCRATCHDIR/$PBS_JOBID, for SGE /tmp/roddyScratch/jobid

NP_FLAGSTATS_IN=${localScratchDirectory}/np_flagstats_in_merged
NP_READBINS_IN=${localScratchDirectory}/np_readbins_in_merged
NP_COVERAGEQC_IN=${localScratchDirectory}/np_coverageqc_in_merged
NP_METRICS_IN=${localScratchDirectory}/np_metrics_in_merged
NP_COMBINEDANALYSIS_IN=${localScratchDirectory}/np_combinedanalysis_in_merged
NP_RAW_BAM_ORG=${localScratchDirectory}/np_raw_bam_original

bamname=`basename ${FILENAME_MERGEDBAM}`

# Ignore chr prefix and set it manually
source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${FILENAME_MERGEDBAM} # Sets CHR_PREFIX and REFERENCE_GENOME

getReadGroupsFromHeader ${FILENAME_MERGEDBAM} # Gets all read groups set in the bam files and saves them in: BAM_READ_GROUPS array

# create empty arrays for per read group named pipes
NP_RAW_BAM_TO_SPLIT_ARRAY=() # will get the cat on the raw bam file to be split
NP_FLAGSTATS_ARRAY=()
NP_READBINS_ARRAY=()
NP_COVERAGEQC_ARRAY=()
NP_METRICS_ARRAY=()
NP_COMBINEDANALYSIS_ARRAY=()

# create empty arrays for file names
FILENAME_DIFFCHROM_MATRIX_ARRAY=()
FILENAME_DIFFCHROM_STATISTICS_ARRAY=()
FILENAME_ISIZES_MATRIX_ARRAY=()
FILENAME_ISIZES_STATISTICS_ARRAY=()
FILENAME_EXTENDED_FLAGSTATS_ARRAY=()
FILENAME_FLAGSTATS_ARRAY=()
FILENAME_READBINS_COVERAGE_ARRAY=()
FILENAME_GENOME_COVERAGE_ARRAY=()
FILENAME_METRICS_ARRAY=()
FILENAME_QCSUMMARY_ARRAY=()
FILENAME_ISIZES_PLOT_ARRAY=()
FILENAME_DIFFCHROM_PLOT_ARRAY=()

for RG in ${BAM_READ_GROUPS[@]}
do
	# Create arrays with named pipes for all read groups
	NP_RAW_BAM_TO_SPLIT_ARRAY[${#NP_RAW_BAM_TO_SPLIT_ARRAY[@]}]=${localScratchDirectory}/raw_bam_${RG}
	NP_FLAGSTATS_ARRAY[${#NP_FLAGSTATS_ARRAY[@]}]=${localScratchDirectory}/np_flagstats_in_${RG}
	NP_READBINS_ARRAY[${#NP_READBINS_ARRAY[@]}]=${localScratchDirectory}/np_readbins_in_${RG}
	NP_COVERAGEQC_ARRAY[${#NP_COVERAGEQC_ARRAY[@]}]=${localScratchDirectory}/np_coverageqc_in_${RG}
	NP_METRICS_ARRAY[${#NP_METRICS_ARRAY[@]}]=${localScratchDirectory}/np_metrics_in_${RG}
	NP_COMBINEDANALYSIS_ARRAY[${#NP_COMBINEDANALYSIS_ARRAY[@]}]=${localScratchDirectory}/np_combinedanalysis_in_${RG}
	# Create arrays with file names
	FILENAME_DIFFCHROM_MATRIX_ARRAY[${#FILENAME_DIFFCHROM_MATRIX_ARRAY[@]}]=`dirname ${FILENAME_DIFFCHROM_MATRIX}`/${SAMPLE}_${RG}.bam_DiffChroms.txt
	FILENAME_DIFFCHROM_STATISTICS_ARRAY[${#FILENAME_DIFFCHROM_STATISTICS_ARRAY[@]}]=`dirname ${FILENAME_DIFFCHROM_STATISTICS}`/${SAMPLE}_${RG}.bam_DiffChroms.png_qcValues.txt
	FILENAME_ISIZES_MATRIX_ARRAY[${#FILENAME_ISIZES_MATRIX_ARRAY[@]}]=`dirname ${FILENAME_ISIZES_MATRIX}`/${SAMPLE}_${RG}_insertsizes.txt
	FILENAME_ISIZES_STATISTICS_ARRAY[${#FILENAME_ISIZES_STATISTICS_ARRAY[@]}]=`dirname ${FILENAME_ISIZES_STATISTICS}`/${SAMPLE}_${RG}_insertsize_plot.png_qcValues.txt
	FILENAME_EXTENDED_FLAGSTATS_ARRAY[${#FILENAME_EXTENDED_FLAGSTATS_ARRAY[@]}]=`dirname ${FILENAME_EXTENDED_FLAGSTATS}`/${SAMPLE}_${RG}.bam_extendedFlagstats.txt
	FILENAME_FLAGSTATS_ARRAY[${#FILENAME_FLAGSTATS_ARRAY[@]}]=`dirname ${FILENAME_FLAGSTATS}`/${SAMPLE}_${RG}.bam_flagstats.txt
	FILENAME_READBINS_COVERAGE_ARRAY[${#FILENAME_READBINS_COVERAGE_ARRAY[@]}]=`dirname ${FILENAME_READBINS_COVERAGE}`/${SAMPLE}_${RG}_readCoverage${FILENAME_READBINS_COVERAGE#*readCoverage}
	FILENAME_GENOME_COVERAGE_ARRAY[${#FILENAME_GENOME_COVERAGE_ARRAY[@]}]=`dirname ${FILENAME_GENOME_COVERAGE}`/${SAMPLE}_${RG}.DepthOfCoverage_Genome.txt
	FILENAME_METRICS_ARRAY[${#FILENAME_METRICS_ARRAY[@]}]=`dirname ${FILENAME_METRICS}`/${SAMPLE}_${RG}.bam.dupmark_metrics.txt
	FILENAME_QCSUMMARY_ARRAY[${#FILENAME_QCSUMMARY_ARRAY[@]}]=`dirname ${FILENAME_QCSUMMARY}`/${SAMPLE}_${RG}_commonBam_wroteQcSummary.txt
	FILENAME_ISIZES_PLOT_ARRAY[${#FILENAME_ISIZES_PLOT_ARRAY[@]}]=`dirname ${FILENAME_ISIZES_PLOT}`/${SAMPLE}_${RG}_insertsize_plot.png
	FILENAME_DIFFCHROM_PLOT_ARRAY[${#FILENAME_DIFFCHROM_PLOT_ARRAY[@]}]=`dirname ${FILENAME_DIFFCHROM_PLOT}`/${SAMPLE}_${RG}.bam_DiffChroms.png
done

mkfifo ${NP_FLAGSTATS_IN} ${NP_READBINS_IN} ${NP_COVERAGEQC_IN} ${NP_METRICS_IN} ${NP_COMBINEDANALYSIS_IN} ${NP_FLAGSTATS_ARRAY[@]} ${NP_READBINS_ARRAY[@]} ${NP_COVERAGEQC_ARRAY[@]} ${NP_METRICS_ARRAY[@]} ${NP_COMBINEDANALYSIS_ARRAY[@]} ${NP_RAW_BAM_TO_SPLIT_ARRAY[@]} ${NP_RAW_BAM_ORG}

MBUF_100M="${MBUFFER_BINARY} -m 100m -q -l /dev/null"
MBUF_2G="${MBUFFER_BINARY} -m 2g -q -l /dev/null"
MBUF_4G="${MBUFFER_BINARY} -m 4g -q -l /dev/null"

# Fill all name pipes with the original bam stream
(cat ${FILENAME_MERGEDBAM} | ${MBUF_2G} | tee ${NP_RAW_BAM_TO_SPLIT_ARRAY[@]} > ${NP_RAW_BAM_ORG}) & procIDcat=$!

# Create tree of input files/named pipes for the merged bam file
(cat ${NP_RAW_BAM_ORG} | tee ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} | ${SAMTOOLS_BINARY} view - | tee ${NP_METRICS_IN} > ${NP_COMBINEDANALYSIS_IN}) & procIDSAMpipe=$!; procIDBBB=$procIDSAMpipe

# Create tree of input files/named pipes for each read group
procIDSAMpipeRG=()
for j in `seq ${#BAM_READ_GROUPS[@]}`
do
	i=`expr $j - 1`
	(cat ${NP_RAW_BAM_TO_SPLIT_ARRAY[$i]} | ${SAMTOOLS_BINARY} view -u -r ${BAM_READ_GROUPS[$i]} - |tee ${NP_FLAGSTATS_ARRAY[$i]} ${NP_COVERAGEQC_ARRAY[$i]} ${NP_READBINS_ARRAY[$i]} | ${SAMTOOLS_BINARY} view - | tee ${NP_METRICS_ARRAY[$i]} > ${NP_COMBINEDANALYSIS_ARRAY[$i]}) & procIDSAMpipeRG[$i]=$!
done

# Create fake flag stats from Barbara for the merged bam file
if [[ -e ${FILENAME_METRICS} ]]
then
    (cat ${NP_METRICS_IN} > /dev/null) & procIDFakeMetrics=$!
else
    (cat ${NP_METRICS_IN} | ${PERL_BINARY} ${TOOL_FAKE_DUPMARK_METRICS} - ${SAMPLE}_${pid} > ${FILENAME_METRICS}) & procIDFakeMetrics=$!
fi

# Create fake flag stats from Barbara for each read group
procIDFakeMetricsRG=()
for j in `seq ${#BAM_READ_GROUPS[@]}`
do
	i=`expr $j - 1`
	if [[ -e ${FILENAME_METRICS_ARRAY[$i]} ]]
	then
		(cat ${NP_METRICS_ARRAY[$i]} > /dev/null) & procIDFakeMetricsRG[$i]=$!
	else
		(cat ${NP_METRICS_ARRAY[$i]} | ${PERL_BINARY} ${TOOL_FAKE_DUPMARK_METRICS} - ${SAMPLE}_${pid} > ${FILENAME_METRICS_ARRAY[$i]}) & procIDFakeMetricsRG[$i]=$!
	fi
done

# in all cases:
# For the merged bam file
# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN} -a ${FILENAME_DIFFCHROM_MATRIX}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX}.tmp  -f ${FILENAME_EXTENDED_FLAGSTATS}.tmp  -m ${FILENAME_ISIZES_STATISTICS}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBA=$!

# make flagstats of piped BAM
(${SAMTOOLS_BINARY} flagstat ${NP_FLAGSTATS_IN} > ${FILENAME_FLAGSTATS}.tmp) & procIDSamtoolsFlagstat=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN} --outputFile=${FILENAME_GENOME_COVERAGE}.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage=$!

# this part often fails with broken pipe, ?? where this comes from. The mbuffer did not help, maybe --processors=4 does?
(${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE} | $MBUF_100M | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}.tmp) & procIDReadbinsCoverage=$!

# Make it per read group...
procIDCBARG=()
procIDSamtoolsFlagstatRG=()
procIDGenomeCoverageRG=()
procIDReadbinsCoverageRG=()
for j in `seq ${#BAM_READ_GROUPS[@]}`
do
	i=`expr $j - 1`
	# SAM output is piped to perl script that calculates various QC measures
	(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_ARRAY[$i]} -a ${FILENAME_DIFFCHROM_MATRIX_ARRAY[$i]}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX_ARRAY[$i]}.tmp  -f ${FILENAME_EXTENDED_FLAGSTATS_ARRAY[$i]}.tmp  -m ${FILENAME_ISIZES_STATISTICS_ARRAY[$i]}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS_ARRAY[$i]}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBARG[$i]=$!

	# make flagstats of piped BAM
	(${SAMTOOLS_BINARY} flagstat ${NP_FLAGSTATS_ARRAY[$i]} > ${FILENAME_FLAGSTATS_ARRAY[$i]}.tmp) & procIDSamtoolsFlagstatRG[$i]=$!

	# genome coverage (depth of coverage and other QC measures in one file)
	(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_ARRAY[$i]} --outputFile=${FILENAME_GENOME_COVERAGE_ARRAY[$i]}.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverageRG[$i]=$!

	# this part often fails with broken pipe, ?? where this comes from. The mbuffer did not help, maybe --processors=4 does?
	(${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_ARRAY[$i]} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE} | $MBUF_100M | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE_ARRAY[$i]}.tmp) & procIDReadbinsCoverageRG[$i]=$!
done

# For the merged bam file
# Some waits for parallel processes. This also depends on the used merge binary.
wait $procIDcat; [[ $? -gt 0 ]] && echo "Error from cat on bamfile" && exit 8
wait $procIDBBB; [[ $? -gt 0 ]] && echo "Error from cat on bamfile" && exit 10
wait $procIDSAMpipe; [[ $? -gt 0 ]] && echo "Error from samtools SAM pipe" && exit 9
wait $procIDReadbinsCoverage; [[ $? -gt 0 ]] && echo "Error from genomeCoverage read bins" && exit 10
wait $procIDGenomeCoverage; [[ $? -gt 0 ]] && echo "Error from coverageQCD" && exit 11
wait $procIDSamtoolsFlagstat; [[ $? -gt 0 ]] && echo "Error from samtools flagstats" && exit 12
wait $procIDCBA; [[ $? -gt 0 ]] && echo "Error from combined QC perl script" && exit 13
wait $procIDFakeMetrics; [[ $? -gt 0 ]] && echo "Error in fake metrics script" && exit 18

mv ${FILENAME_DIFFCHROM_MATRIX}.tmp ${FILENAME_DIFFCHROM_MATRIX}
mv ${FILENAME_ISIZES_MATRIX}.tmp ${FILENAME_ISIZES_MATRIX}
mv ${FILENAME_EXTENDED_FLAGSTATS}.tmp ${FILENAME_EXTENDED_FLAGSTATS}
mv ${FILENAME_ISIZES_STATISTICS}.tmp ${FILENAME_ISIZES_STATISTICS}
mv ${FILENAME_DIFFCHROM_STATISTICS}.tmp ${FILENAME_DIFFCHROM_STATISTICS}
mv ${FILENAME_FLAGSTATS}.tmp ${FILENAME_FLAGSTATS}
mv ${FILENAME_READBINS_COVERAGE}.tmp ${FILENAME_READBINS_COVERAGE}
mv ${FILENAME_GENOME_COVERAGE}.tmp ${FILENAME_GENOME_COVERAGE}

# For each read group
for j in `seq ${#BAM_READ_GROUPS[@]}`
do
	i=`expr $j - 1`
	wait ${procIDSAMpipeRG[$i]}; [[ $? -gt 0 ]] && echo "Error from samtools SAM pipe ${BAM_READ_GROUPS[$i]}" && exit 19
	wait ${procIDReadbinsCoverageRG[$i]}; [[ $? -gt 0 ]] && echo "Error from genomeCoverage read bins ${BAM_READ_GROUPS[$i]}" && exit 20
	wait ${procIDGenomeCoverageRG[$i]}; [[ $? -gt 0 ]] && echo "Error from coverageQCD ${BAM_READ_GROUPS[$i]}" && exit 21
	wait ${procIDSamtoolsFlagstatRG[$i]}; [[ $? -gt 0 ]] && echo "Error from samtools flagstats ${BAM_READ_GROUPS[$i]}" && exit 22
	wait ${procIDCBARG[$i]}; [[ $? -gt 0 ]] && echo "Error from combined QC perl script ${BAM_READ_GROUPS[$i]}" && exit 23
	wait ${procIDFakeMetricsRG[$i]}; [[ $? -gt 0 ]] && echo "Error in fake metrics script ${BAM_READ_GROUPS[$i]}" && exit 24

	mv ${FILENAME_DIFFCHROM_MATRIX_ARRAY[$i]}.tmp ${FILENAME_DIFFCHROM_MATRIX_ARRAY[$i]}
	mv ${FILENAME_ISIZES_MATRIX_ARRAY[$i]}.tmp ${FILENAME_ISIZES_MATRIX_ARRAY[$i]}
	mv ${FILENAME_EXTENDED_FLAGSTATS_ARRAY[$i]}.tmp ${FILENAME_EXTENDED_FLAGSTATS_ARRAY[$i]}
	mv ${FILENAME_ISIZES_STATISTICS_ARRAY[$i]}.tmp ${FILENAME_ISIZES_STATISTICS_ARRAY[$i]}
	mv ${FILENAME_DIFFCHROM_STATISTICS_ARRAY[$i]}.tmp ${FILENAME_DIFFCHROM_STATISTICS_ARRAY[$i]}
	mv ${FILENAME_FLAGSTATS_ARRAY[$i]}.tmp ${FILENAME_FLAGSTATS_ARRAY[$i]}
	mv ${FILENAME_READBINS_COVERAGE_ARRAY[$i]}.tmp ${FILENAME_READBINS_COVERAGE_ARRAY[$i]}
	mv ${FILENAME_GENOME_COVERAGE_ARRAY[$i]}.tmp ${FILENAME_GENOME_COVERAGE_ARRAY[$i]}
done

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

# For the merged bam file
[[ -f ${FILENAME_QCSUMMARY}_WARNINGS.txt ]] && rm ${FILENAME_QCSUMMARY}_WARNINGS.txt

${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY -p $PID -s $SAMPLE -r all_merged -l $analysis_type -w ${FILENAME_QCSUMMARY}_WARNINGS.txt -f $FILENAME_FLAGSTATS -d $FILENAME_DIFFCHROM_STATISTICS -i $FILENAME_ISIZES_STATISTICS -c $FILENAME_GENOME_COVERAGE -m ${FILENAME_METRICS} > ${FILENAME_QCSUMMARY}_temp && mv ${FILENAME_QCSUMMARY}_temp $FILENAME_QCSUMMARY || ( echo "Error from writeQCsummary.pl" && exit 14)


# For each read group
for j in `seq ${#BAM_READ_GROUPS[@]}`
do
	i=`expr $j - 1`
	[[ -f ${FILENAME_QCSUMMARY_ARRAY[$i]}_WARNINGS.txt ]] && rm ${FILENAME_QCSUMMARY_ARRAY[$i]}_WARNINGS.txt

	${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY -p $PID -s $SAMPLE -r all_merged -l $analysis_type -w ${FILENAME_QCSUMMARY_ARRAY[$i]}_WARNINGS.txt -f ${FILENAME_FLAGSTATS_ARRAY[$i]} -d ${FILENAME_DIFFCHROM_STATISTICS_ARRAY[$i]} -i ${FILENAME_ISIZES_STATISTICS_ARRAY[$i]} -c ${FILENAME_GENOME_COVERAGE_ARRAY[$i]} -m ${FILENAME_METRICS_ARRAY[$i]} > ${FILENAME_QCSUMMARY_ARRAY[$i]}_temp && mv ${FILENAME_QCSUMMARY_ARRAY[$i]}_temp ${FILENAME_QCSUMMARY_ARRAY[$i]} || ( echo "Error from writeQCsummary.pl" && exit 14)

done

[[ ${useSingleEndProcessing-false} == true ]] && exit 0

# For the merged bam file
${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAME_ISIZES_MATRIX} ${FILENAME_ISIZES_STATISTICS} ${FILENAME_ISIZES_PLOT}_temp "PE insertsize of ${bamname} (rmdup)" && mv  ${FILENAME_ISIZES_PLOT}_temp ${FILENAME_ISIZES_PLOT} || ( echo "Error from insert sizes plotter" && exit 22)
${RSCRIPT_BINARY} ${TOOL_PLOT_DIFFCHROM} -i "$FILENAME_DIFFCHROM_MATRIX" -s "$FILENAME_DIFFCHROM_STATISTICS" -o "${FILENAME_DIFFCHROM_PLOT}_temp" && mv  ${FILENAME_DIFFCHROM_PLOT}_temp ${FILENAME_DIFFCHROM_PLOT} || ( echo "Error from chrom_diff.r" && exit 23)

# For each read group
for j in `seq ${#BAM_READ_GROUPS[@]}`
do
	i=`expr $j - 1`
	${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAME_ISIZES_MATRIX_ARRAY[$i]} ${FILENAME_ISIZES_STATISTICS_ARRAY[$i]} ${FILENAME_ISIZES_PLOT_ARRAY[$i]}_temp "PE insertsize of ${BAM_READ_GROUPS[$i]}" && mv  ${FILENAME_ISIZES_PLOT_ARRAY[$i]}_temp ${FILENAME_ISIZES_PLOT_ARRAY[$i]} || ( echo "Error from insert sizes plotter ${BAM_READ_GROUPS[$i]}" && exit 22)

	${RSCRIPT_BINARY} ${TOOL_PLOT_DIFFCHROM} -i "${FILENAME_DIFFCHROM_MATRIX_ARRAY[$i]}" -s "${FILENAME_DIFFCHROM_STATISTICS_ARRAY[$i]}" -o "${FILENAME_DIFFCHROM_PLOT_ARRAY[$i]}_temp" && mv  ${FILENAME_DIFFCHROM_PLOT_ARRAY[$i]}_temp ${FILENAME_DIFFCHROM_PLOT_ARRAY[$i]} || ( echo "Error from chrom_diff.r ${BAM_READ_GROUPS[$i]}" && exit 23)
done
