#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

#PBS -l nodes=1:ppn=1
#PBS -l mem=750m
#PBS -l walltime=00:10:00
#PBS -m a

## the outfile name structure is necessary for the R script: sample_pid_readCoverage_${WINDOW_SIZE}kb_windows.txt
#OUTFILENAME=$DIR_COVERAGE/tumor_AvD_59208_readCoverage_${WINDOW_SIZE}kb_windows.txt

# plotting: eigener Job, nodes=1:ppn=1,walltime=00:20:00,mem=1g
# control vs. tumor (dependency auf alleset -xv genomeCoverage.py-Jobs)
# mehrere: am besten in einem Job: control vs. tumor_01 control vs. tumor_02, control vs. tumor_03, ... (control immer als erstes Argument)

#$DIR_COVERAGE/control_AvD_59208_vs_tumor_AvD_59208_readCoverage_${WINDOW_SIZE}kb_windows.png

# plot genome-wide first sample (control) against second sample (tumor)
# If FILENAME_TUMOR is set then we try to plot for two samples
if [[ "x" == ${FILENAME_TUMOR+x} ]]; then
    ${RSCRIPT_BINARY} ${TOOL_COVERAGE_PLOT_SCRIPT} "${FILENAME_CONTROL}" "${FILENAME_TUMOR}" "read" "${WINDOW_SIZE}" "${PLOTFILE}"
else
    ${RSCRIPT_BINARY} ${TOOL_COVERAGE_PLOT_SCRIPT} "${FILENAME_CONTROL}" "read" "${WINDOW_SIZE}" "${PLOTFILE}"
fi
