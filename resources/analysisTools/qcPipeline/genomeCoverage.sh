#!/bin/bash

#PBS -l nodes=1:ppn=4
#PBS -l walltime=5:00:00
#PBS -m n

source ${CONFIG_FILE}

attachement=""
if [[ ${FILENAME} == *targetExtract* ]] && [[ -n ${TARGETSIZE} ]]
then
	attachement="--targetsize=${TARGETSIZE}"
fi

${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${FILENAME} --outputFile=${FILENAME_COVERAGE}_tmp --processors=4 --basequalCutoff=$BASE_QUALITY_CUTOFF --ungappedSizes=${CHROM_SIZES_FILE} $attachement && mv ${FILENAME_COVERAGE}_tmp ${FILENAME_COVERAGE} || (echo "Non-zero exit code from coverageQc" && exit 2)

