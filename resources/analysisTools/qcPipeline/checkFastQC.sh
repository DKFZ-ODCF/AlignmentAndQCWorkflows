#!/bin/bash

#PBS -l walltime=10:00:00
#PBS -l nodes=1
#PBS -m a

set -uvex

source "$TOOL_BASH_LIB"
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

    mv $TMP_DIR/*.zip $FILENAME_FASTQC
    rm -r "$TMP_DIR"
fi
