#!/bin/bash

#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=2
#PBS -m a
#PBS -l mem=8g
#PBS -j oe


# no BAM file is produced, everything is piped
# there are no readbins for the target extracted regions
# the parent BAM's coverageQC file is needed

source ${CONFIG_FILE}

set -o pipefail

# TODO: make this come from outside
TARGETS_PLOT=`dirname ${FILENAME_TARGETS_WITH_COVERAGE_TEXT}`"/${SAMPLE}_${PID}_targetCovDistribution.png"

MBUF_1G="${MBUFFER_BINARY} -m 1g -q -l /dev/null"

localScratchDirectory=${PBS_SCRATCH_DIR}/${PBS_JOBID}
NP_COMP_IN=${localScratchDirectory}/np_compression_in
NP_FLAGSTATS_IN=${localScratchDirectory}/np_flagstats_in
NP_COVERAGEQC_IN=${localScratchDirectory}/np_coverageqc_in
NP_COMBINEDANALYSIS_IN=${localScratchDirectory}/np_combinedanalysis_in

mkfifo ${NP_COMP_IN} ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} ${NP_COMBINEDANALYSIS_IN}

samplepid=${SAMPLE}_${PID}

# make a SAM pipe for the Perl tool
${SAMTOOLS_BINARY} view ${NP_COMP_IN} | ${MBUF_1G} > ${NP_COMBINEDANALYSIS_IN} & procIDSAMpipe=$!

# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN} -a ${FILENAME_DIFFCHROM_MATRIX}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX}.tmp -f ${FILENAME_EXTENDED_FLAGSTATS}.tmp  -m ${FILENAME_ISIZES_STATISTICS}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBA=$!

# make flagstats of piped BAM
(${SAMTOOLS_BINARY} flagstat ${NP_FLAGSTATS_IN} > ${FILENAME_FLAGSTATS}.tmp) & procIDSamtoolsFlagstat=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN} --outputFile=${FILENAME_GENOME_COVERAGE}.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE} --targetsize=${TARGETSIZE}) & procIDGenomeCoverage=$!

# create a streamed BAM file of target regions with intersectBed
# on this, also calculate the per-target coverage with coverageBed per base and parse per perl script per region
(set -o pipefail; $INTERSECTBED_BINARY $INTERSECTBED_OPTIONS -abam $FILENAME_PARENTBAM -b $TARGET_REGIONS_FILE | ${MBUF_1G} | \
tee ${NP_COMP_IN} ${NP_FLAGSTATS_IN} ${NP_COVERAGEQC_IN} | \
$SAMTOOLS_BINARY view $SAMTOOLS_VIEW_OPTIONS - | ${MBUF_1G} | \
$COVERAGEBED_BINARY $COVERAGEBED_OPTIONS -abam stdin -b $TARGET_REGIONS_FILE | \
${PERL_BINARY} ${TOOL_TARGET_COVERAGE_PERL_SCRIPT} - > ${FILENAME_TARGETS_WITH_COVERAGE_TEXT}.tmp ; \
echo $? > ${DIR_TEMP}/${samplepid}_ec_target) & procIDtargetExtr=$!

wait $procIDtargetExtr; [[ ! `cat ${DIR_TEMP}/${samplepid}_ec_target` -eq "0" ]] && echo "intersectBed - samtools pipe returned a non-zero exit code and the job will die now." && exit 100

wait $procIDGenomeCoverage; [[ $? -gt 0 ]] && echo "Error from coverageQCD" && exit 11
wait $procIDSamtoolsFlagstat; [[ $? -gt 0 ]] && echo "Error from samtools flagstats" && exit 12
wait $procIDCBA; [[ $? -gt 0 ]] && echo "Error from combined QC perl script" && exit 13


mv ${FILENAME_DIFFCHROM_MATRIX}.tmp ${FILENAME_DIFFCHROM_MATRIX}
mv ${FILENAME_ISIZES_MATRIX}.tmp ${FILENAME_ISIZES_MATRIX}
mv ${FILENAME_EXTENDED_FLAGSTATS}.tmp ${FILENAME_EXTENDED_FLAGSTATS}
mv ${FILENAME_ISIZES_STATISTICS}.tmp ${FILENAME_ISIZES_STATISTICS}
mv ${FILENAME_DIFFCHROM_STATISTICS}.tmp ${FILENAME_DIFFCHROM_STATISTICS}
mv ${FILENAME_GENOME_COVERAGE}.tmp ${FILENAME_GENOME_COVERAGE}
mv ${FILENAME_FLAGSTATS}.tmp ${FILENAME_FLAGSTATS}
mv ${FILENAME_TARGETS_WITH_COVERAGE_TEXT}.tmp ${FILENAME_TARGETS_WITH_COVERAGE_TEXT}


# QC summary gets the parent coverageQC file $FILENAME_PARENTBAM_COVERAGE
# and the target extract QC files created above
# therefore the parameters are slightly different (-c and -t)
# FILENAME_GENOME_COVERAGE will be the coverageQC of the target extracted BAM, with on target coverage in the last line
(${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY -p $PID -s $SAMPLE -r all_merged -l exome_target -w ${FILENAME_QCSUMMARY}_WARNINGS.txt -f $FILENAME_FLAGSTATS -d $FILENAME_DIFFCHROM_STATISTICS -i $FILENAME_ISIZES_STATISTICS -c $FILENAME_PARENTBAM_COVERAGE -t $FILENAME_GENOME_COVERAGE > ${FILENAME_QCSUMMARY}.tmp && mv ${FILENAME_QCSUMMARY}.tmp $FILENAME_QCSUMMARY) || ( echo "Error from writeQCsummary.pl" && exit 14)

# plot the cumulative target coverage
${RSCRIPT_BINARY} ${TOOL_ON_TARGET_COVERAGE_PLOTTER_BINARY} ${FILENAME_TARGETS_WITH_COVERAGE_TEXT} ${TARGETS_PLOT}.tmp "${samplepid}" && mv ${TARGETS_PLOT}.tmp ${TARGETS_PLOT} || ( echo "Error from on target coverage plotter" && exit 15)

# if the BAM only contains single end reads, there can be no pairs to have insert sizes or ends mapping to different chromosomes
# and plotting would fail
grep -w "NA" $FILENAME_DIFFCHROM_STATISTICS && useSingleEndProcessing=true

[[ ${useSingleEndProcessing-false} == true ]] && exit 0

# plot the insert size distribution
${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAME_ISIZES_MATRIX} ${FILENAME_ISIZES_STATISTICS} ${FILENAME_ISIZES_PLOT}_temp "PE insertsize of ${FILEINFO} (rmdup)" && mv  ${FILENAME_ISIZES_PLOT}_temp ${FILENAME_ISIZES_PLOT} || ( echo "Error from insert sizes plotter" && exit 16)

# plot the paired end aberrations (reads on different chromosomes)
${RSCRIPT_BINARY} ${TOOL_PLOT_DIFFCHROM} -i "$FILENAME_DIFFCHROM_MATRIX" -s "$FILENAME_DIFFCHROM_STATISTICS" -o "${FILENAME_DIFFCHROM_PLOT}.tmp" && mv ${FILENAME_DIFFCHROM_PLOT}.tmp ${FILENAME_DIFFCHROM_PLOT} || ( echo "Error from chromdiff plotter" && exit 17)
