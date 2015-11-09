#!/bin/bash

#PBS -l nodes=1:ppn=2
#PBS -l walltime=2:00:00
#PBS -m a
#PBS -l mem=4g
#PBS -j oe

source ${CONFIG_FILE}

readlen=`$SAMTOOLS_BINARY view ${FILENAME} | head -n 1 | awk '{print length($10)}'`

$SAMTOOLS_BINARY view $SAMTOOLS_VIEW_OPTIONS ${FILENAME} | $COVERAGEBED_BINARY $COVERAGEBED_OPTIONS -abam stdin -b ${TARGET_REGIONS_FILE} | perl -p -e 'BEGIN {$rl=shift} chomp; @a=split(/\t/); $cov=$a[-4]*$rl/($a[2]-$a[1]); $_ .= "\t$cov\n" ' $readlen > $TARGETS_WITH_COVERAGE_TEXT

# Is now in createOnTargetCoveragePlot.sh R -f ${ONTARGETCOVERAGEPLOTTER} --no-save --no-restore --args $TARGETS_WITH_COVERAGE_TEXT $TARGETS_PLOT "${FILENAME_PREFIX}"
