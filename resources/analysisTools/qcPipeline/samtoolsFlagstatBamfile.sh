#!/bin/bash

#PBS -l walltime=5:00:00
#PBS -l nodes=1
#PBS -l mem=50m
#PBS -m a

#[ $DEBUG = 'TRUE' ] && set -x && set -v

TMP_FILE=${FILENAME_FLAGSTAT}_temp

${SAMTOOLS_BINARY} flagstat ${FILENAME} > ${TMP_FILE}

mv ${TMP_FILE} ${FILENAME_FLAGSTAT}
