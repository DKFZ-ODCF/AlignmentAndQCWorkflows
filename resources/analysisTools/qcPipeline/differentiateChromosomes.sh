#!/bin/bash

#PBS -l walltime=5:00:00
#PBS -l nodes=1
#PBS -l mem=1g
#PBS -m a

filenameTemp=${FILENAMED}_tmp
${PYTHON_BINARY} ${TOOL_PAIRED_END_READ_ABERRATIONS_SCRIPT} --alignmentFile=${FILENAME} --outputFile=${filenameTemp} --printSummaryToStdOut --plotDistribution --Rpath=`dirname ${TOOL_PAIRED_END_READ_ABERRATIONS_SCRIPT}`

[[ $? == 0 ]] && mv ${filenameTemp} ${FILENAMED} && exit 0

echo "Error in script"
exit 5