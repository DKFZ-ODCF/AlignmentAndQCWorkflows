#!/bin/bash

#PBS -l walltime=10:00:00
#PBS -l nodes=1
#PBS -m a

source ${CONFIG_FILE}

# Check for single lane processing.
# For single lanes an empty fastqc file is created
[[ -z "${useSingleEndProcessing+x}" && ! -f ${RAW_SEQ} ]] && touch ${FILENAME_FASTQC} && exit 0

# by default, fastqc writes a zipped file with the infile name as prefix to the outdir

DIRNAME_FASTQC=`dirname $FILENAME_FASTQC`

TMP_DIR=$DIRNAME_FASTQC/$PBS_JOBID
mkdir $TMP_DIR

${FASTQC_BINARY} ${RAW_SEQ} --noextract -o $TMP_DIR

mv $TMP_DIR/* $FILENAME_FASTQC

rm -r $TMP_DIR
