#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

#PBS -l walltime=5:00:00
#PBS -l nodes=1
#PBS -l mem=300m
#PBS -m a

set -o pipefail
FILENAMESHORT=`basename ${FILENAME}`
FILEINFO=${FILENAMESHORT%_*}


#${SAMTOOLS_BINARY} view -f 67 -F 4 ${FILENAME} | awk '{ printf "%s\n",$9 }'> ${FILENAMED}
#R -f $TOOLSDIR/insertSizeDistribution.r --no-save --no-restore --args "${FILENAMED}" "${FILENAMEP}" "PE insertsize of ${FILEINFO} (rmdup)"
# pipe insert sizes of proper pairs to R for getting median, standard deviation and bimodal plot
#${SAMTOOLS_BINARY} view -f 67 -F 4 ${FILENAME} | awk '{ printf "%s\n",$9 }' | Rscript $TOOLSDIR/insertSizeDistribution.r "${FILENAMED}" "${FILENAMEP}" "PE insertsize of ${FILEINFO} (rmdup)"

# -F 1024 because duplicates have to be disregarded if they are only marked
# -F 4 is for excluding unmapped reads, but these are not proper pair (-f 67) anyways!

# bwa mem seems to give flag "proper pair" to each read => insert size can go over 200 million bp, makes no sense at all
# Thus take the upper limit from BWA_SAMPESORT_OPTIONS="-a 1000" or default 1000

isizelimit=1000
if [[ "x" == ${BWA_SAMPESORT_OPTIONS+x} ]]
then
	isizelimit=`echo ${BWA_SAMPESORT_OPTIONS} | cut -d " " -f 2`
fi

${SAMTOOLS_BINARY} view -f 67 -F 1024 ${FILENAME} | awk '{ printf "%s\n",$9 }' | ${PERL_BINARY} ${TOOL_INSERT_SIZES_BUCKET_SORT_SCRIPT} - ${FILENAMEP}_qcValues.txt.tmp $isizelimit > ${FILENAMED}.tmp

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the samtools view-perl pipe"
	exit 21
else
	mv ${FILENAMEP}_qcValues.txt.tmp ${FILENAMEP}_qcValues.txt
	mv ${FILENAMED}.tmp ${FILENAMED}
fi

${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAMED} ${FILENAMEP}_qcValues.txt ${FILENAMEP} "PE insertsize of ${FILEINFO} (rmdup)"


