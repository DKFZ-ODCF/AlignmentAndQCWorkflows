#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

#PBS -l nodes=1:ppn=8
#PBS -l mem=1500m
#PBS -l walltime=12:00:00
#PBS -m a

set -o pipefail

#WINDOW_SIZE=1
#TOOLSDIR=/home/hutter/newCoverage
#DIR_COVERAGE=/icgc/dkfzlsdf/analysis/hutter/AvD/results_per_pid/AvD_59208/coverage_new
#DIR_ALIGNMENT=/icgc/dkfzlsdf/analysis/hutter/AvD/results_per_pid/AvD_59208/alignment
#FILENAME=$DIR_ALIGNMENT/tumor_AvD_59208_merged.bam.rmdup.bam
## the outfile name structure is necessary for the R script: sample_pid_readCoverage_${WINDOW_SIZE}kb_windows.txt
#OUTFILENAME=$DIR_COVERAGE/tumor_AvD_59208_readCoverage_${WINDOW_SIZE}kb_windows.txt

TMP_FILE=${OUTFILENAME}_tmp

# Python script counts reads for all contigs
[[ ${useDefaultCoverageQCTool-false} == true ]] &&  ${PYTHON_BINARY} ${TOOL_GENOME_COVERAGE_PY_SCRIPT} --alignmentFile=${FILENAME},${TMP_FILE} --processors=12 --countReads --windowSize=${WINDOW_SIZE} --ignore_chrRandom_chrM_hap_chrUn
[[ ${useDefaultCoverageQCTool-false} == false ]] && ${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${FILENAME} --outputFile=/dev/stdout --processors=2 --mode=countReads --windowSize=${WINDOW_SIZE} | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${TMP_FILE}

[[ $? == 0 ]] && mv ${TMP_FILE} ${OUTFILENAME} && exit 0

echo "Error in readbins genome coverage"
exit 5
