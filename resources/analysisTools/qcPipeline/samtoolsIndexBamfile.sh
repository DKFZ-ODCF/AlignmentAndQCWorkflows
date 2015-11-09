#!/bin/bash

#PBS -l walltime=5:00:00
#PBS -l nodes=1
#PBS -l mem=50m
#PBS -m a

source ${CONFIG_FILE}

TMP_FILE=${FILENAME}_TMP
LOG_FILE=${FILENAME}_ERRLOG
IDX_FILE=${FILENAME}.bai

${SAMTOOLS_BINARY} index ${FILENAME} ${TMP_FILE} 2> $LOG_FILE
[ -s $LOG_FILE ] && cat $LOG_FILE && return -1

mv ${TMP_FILE} ${IDX_FILE}
