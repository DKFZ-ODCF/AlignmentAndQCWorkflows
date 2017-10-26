#!/bin/bash

#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=2
#PBS -m a
#PBS -l mem=5g
#PBS -j oe

set -o pipefail

LCLSCRATCH=${RODDY_SCRATCH}
TMP_PIPE=${LCLSCRATCH}/NP_IBED
TMP_FILE=${targetOnlyBam}.tmp
TMP_TARGETSCOV=${targetsWithCoverageText}.tmp
mkfifo ${TMP_PIPE}

${SAMTOOLS_BINARY} index ${TMP_PIPE} ${TMP_FILE}.bai &
$INTERSECTBED_BINARY $INTERSECTBED_OPTIONS -abam $rawBam -b $TARGET_REGIONS_FILE | tee $TMP_PIPE $TMP_FILE | $SAMTOOLS_BINARY view $SAMTOOLS_VIEW_OPTIONS - | $COVERAGEBED_BINARY $COVERAGEBED_OPTIONS -abam stdin -b $TARGET_REGIONS_FILE | ${PERL_BINARY} ${TOOL_TARGET_COVERAGE_PERL_SCRIPT} - > ${TMP_TARGETSCOV} && mv ${TMP_TARGETSCOV} ${targetsWithCoverageText}

if [[ "$?" != 0 ]]
then
    echo "intersectBed-samtools view-coverageBed-targetCov pipe returned non-zero exit value; exiting..."
    exit 2
fi

wait

mv ${TMP_FILE} ${targetOnlyBam}
mv ${TMP_FILE}.bai ${targetOnlyBam}.bai
touch ${targetOnlyBam}.bai
