#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

#PBS -l nodes=1:ppn=1
#PBS -l walltime=5:00:00
#PBS -l mem=3g
#PBS -m a

[[ ! -f $DIR_METRICS ]] && mkdir -p $DIR_METRICS

cd $DIR_METRICS

picard.sh CollectMultipleMetrics INPUT=${FILENAME} REFERENCE_SEQUENCE=${INDEX_PREFIX} ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT OUTPUT=${COLLECT_METRICS_PREFIX} ${COLLECT_METRICS_TOOLS}
