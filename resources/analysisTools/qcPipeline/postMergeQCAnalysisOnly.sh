#!/bin/bash

source ${CONFIG_FILE}

set -o pipefail
set -u
today=`date +'%Y-%m-%d_%Hh%M'`
workDirectory=${DIR_TEMP}/${RODDY_JOBID}
localScratchDirectory=${RODDY_SCRATCH}	# is for PBS $PBS_SCRATCHDIR/$PBS_JOBID, for SGE /tmp/roddyScratch/jobid
tempDirectory=${FILENAME_MERGEDBAM}_MOMDUP
mkdir -p $workDirectory

tempBamFile=${tempDirectory}/FILENAME_MERGEDBAM_TMP
tempFlagstatsFile=${FILENAME_FLAGSTATS}.tmp
tempIndexFile=${tempBamFile}.bai.tmp

NP_PIC_OUT=${localScratchDirectory}/np_picard_out.sam
NP_SAM_IN=${localScratchDirectory}/np_compression_in
NP_INDEX_IN=${localScratchDirectory}/np_index_in
NP_FLAGSTATS_IN=${localScratchDirectory}/np_flagstats_in
NP_READBINS_IN=${localScratchDirectory}/np_readbins_in
NP_COVERAGEQC_IN=${localScratchDirectory}/np_coverageqc_in
NP_METRICS_IN=${localScratchDirectory}/np_metrics_in
NP_COMBINEDANALYSIS_IN=${localScratchDirectory}/np_combinedanalysis_in

returnCodeMarkDuplicatesFile=${tempDirectory}/rcMarkDup.txt
bamname=`basename ${FILENAME_MERGEDBAM}`

# Convert Roddy input array to colon separated array.
strlen=`expr ${#INPUT_FILES} - 2`
tempInputFiles=""
for inFile in ${INPUT_FILES:1:$strlen}; do tempInputFiles=${tempInputFiles}":"${inFile}; done
INPUT_FILES=${tempInputFiles:1}

bamFileExists=true

# Ignore chr prefix and set it manually
source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${FILENAME_MERGEDBAM} # Sets CHR_PREFIX and REFERENCE_GENOME

mkfifo ${NP_PIC_OUT} ${NP_SAM_IN} ${NP_INDEX_IN} ${NP_FLAGSTATS_IN} ${NP_READBINS_IN} ${NP_COVERAGEQC_IN} ${NP_METRICS_IN} ${NP_COMBINEDANALYSIS_IN}

MBUF_100M="${MBUFFER_BINARY} -m 100m -q -l /dev/null"
MBUF_2G="${MBUFFER_BINARY} -m 2g -q -l /dev/null"
MBUF_4G="${MBUFFER_BINARY} -m 4g -q -l /dev/null"

(cat ${FILENAME_MERGEDBAM} | ${MBUF_2G} | tee ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} | ${SAMTOOLS_BINARY} view - | tee ${NP_METRICS_IN} > ${NP_COMBINEDANALYSIS_IN}) & procIDSAMpipe=$!; procIDBBB=$procIDSAMpipe

if [[ -e ${FILENAME_METRICS} ]]
then
    (cat ${NP_METRICS_IN} > /dev/null) & procIDFakeMetrics=$!
else
    (cat ${NP_METRICS_IN} | ${PERL_BINARY} ${TOOL_FAKE_DUPMARK_METRICS} - ${SAMPLE}_${pid} > ${FILENAME_METRICS}) & procIDFakeMetrics=$!
fi

# in all cases:
# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN} -a ${FILENAME_DIFFCHROM_MATRIX}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX}.tmp  -f ${FILENAME_EXTENDED_FLAGSTATS}.tmp  -m ${FILENAME_ISIZES_STATISTICS}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBA=$!

# make flagstats of piped BAM
(${SAMTOOLS_BINARY} flagstat ${NP_FLAGSTATS_IN} > ${tempFlagstatsFile}) & procIDSamtoolsFlagstat=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN} --outputFile=${FILENAME_GENOME_COVERAGE}.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage=$!

# this part often fails with broken pipe, ?? where this comes from. The mbuffer did not help, maybe --processors=4 does?
(${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE} | $MBUF_100M | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}.tmp) & procIDReadbinsCoverage=$!

# Some waits for parallel processes. This also depends on the used merge binary.
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
mv ${tempFlagstatsFile} ${FILENAME_FLAGSTATS}
mv ${FILENAME_READBINS_COVERAGE}.tmp ${FILENAME_READBINS_COVERAGE}
mv ${FILENAME_GENOME_COVERAGE}.tmp ${FILENAME_GENOME_COVERAGE}


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

${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY -p $PID -s $SAMPLE -r all_merged -l $analysis_type -w ${FILENAME_QCSUMMARY}_WARNINGS.txt -f $FILENAME_FLAGSTATS -d $FILENAME_DIFFCHROM_STATISTICS -i $FILENAME_ISIZES_STATISTICS -c $FILENAME_GENOME_COVERAGE -m ${FILENAME_METRICS} > ${FILENAME_QCSUMMARY}_temp && mv ${FILENAME_QCSUMMARY}_temp $FILENAME_QCSUMMARY || ( echo "Error from writeQCsummary.pl" && exit 14)

[[ -d $tempDirectory ]] && rm -rf $tempDirectory

[[ ${useSingleEndProcessing-false} == true ]] && exit 0

${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAME_ISIZES_MATRIX} ${FILENAME_ISIZES_STATISTICS} ${FILENAME_ISIZES_PLOT}_temp "PE insertsize of ${bamname} (rmdup)" && mv  ${FILENAME_ISIZES_PLOT}_temp ${FILENAME_ISIZES_PLOT} || ( echo "Error from insert sizes plotter" && exit 22)
${RSCRIPT_BINARY} ${TOOL_PLOT_DIFFCHROM} -i "$FILENAME_DIFFCHROM_MATRIX" -s "$FILENAME_DIFFCHROM_STATISTICS" -o "${FILENAME_DIFFCHROM_PLOT}_temp" && mv  ${FILENAME_DIFFCHROM_PLOT}_temp ${FILENAME_DIFFCHROM_PLOT} || ( echo "Error from chrom_diff.r" && exit 23)

