#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

#PBS -l walltime=10:00:00
#PBS -l nodes=1
#PBS -m a

set -uvex

source "$TOOL_WORKFLOW_LIB"
printInfo

# Check for single lane processing.
# For single lanes an empty fastqc file is created
if [[ -v useSingleEndProcessing && "$useSingleEndProcessing" == true && ! -f ${RAW_SEQ} ]]; then
    touch ${FILENAME_FASTQC}
else
    # by default, fastqc writes a zipped file with the infile name as prefix to the outdir

    DIRNAME_FASTQC=`dirname $FILENAME_FASTQC`
    TMP_DIR=$DIRNAME_FASTQC/$RODDY_JOBID
    mkdir "$TMP_DIR"

    ${FASTQC_BINARY} ${RAW_SEQ} --noextract -o $TMP_DIR \
        || throw 1 "Error during FASTQC"

    "${TOOL_FASTQC_CLASSIFY}" \
        <(unzip -p "$TMP_DIR/*.zip" "*/fastqc_data.txt") \
        > "${FILENAME_FASTQ_QC_STATUS}" || throw 10 "Error classifying the FASTQ quality"

    mv $TMP_DIR/*.zip $FILENAME_FASTQC
    rm -r "$TMP_DIR"
fi
