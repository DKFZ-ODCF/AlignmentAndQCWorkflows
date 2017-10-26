#!/bin/bash

#PBS -l nodes=1:ppn=2
#PBS -l walltime=2:00:00
#PBS -m a
#PBS -l mem=4g
#PBS -j oe

R -f ${TOOL_ON_TARGET_COVERAGE_PLOTTER_BINARY} --no-save --no-restore --args ${TARGETS_WITH_COVERAGE_TEXT} ${TARGETS_PLOT} "${FILENAME_PREFIX}"
