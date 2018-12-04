#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

source $TOOL_BASH_LIB

set -xuv

bamfile=${FILENAME_MERGED_BAM}

# Get atomic filename of bam_file
bamfile_basename=`basename ${bamfile}`;
bamfile_atomic_name="${bamfile_basename%.*}";

meth_call_files_basename=${bamfile_basename};
meth_call_files_dirname=`dirname ${FILENAME_METH_CALLS_META_CHECKPOINT}`;
meth_call_metric_files_dirname=`dirname ${FILENAME_METH_CALLS_META_METRICS_CHECKPOINT}`;
FILENAME_PREFIX_METHYLATION_CALLING=${meth_call_files_dirname}/${meth_call_files_basename};
FILENAME_PREFIX_METHYLATION_CALLING_METRICS=${meth_call_metric_files_dirname}/${meth_call_files_basename}

# Create named pipes in ${RODDY_SCRATCH}
SORTED_CHROMOSOMES=${RODDY_SCRATCH}/sorted_chromosomes.txt;
declare -a CHROMOSOME_INDICES="$CHROMOSOME_INDICES"

# Sort Chromosome Indices by Chromosome length and store them in CHROMOSOME_INDICES_SORTED
# Create associative array containing (chromosome: chromosome_size) pairs
declare -A chromosome_sizes;
while read line; do
	chrom=`echo ${line} | cut -f 1 -d " "`;
	size=`echo ${line} | cut -f 2 -d " "`;
	chromosome_sizes[${chrom}]=${size};
done < ${CHROM_SIZES_FILE};

for chrom in ${CHROMOSOME_INDICES[@]}; do
        echo -e ${chrom}"\t"${chromosome_sizes[${chrom}]};
done | sort -n -k2 >> ${SORTED_CHROMOSOMES};

CHROMOSOME_INDICES_SORTED=();
while read line; do
        elem=`echo ${line} | cut -f 1 -d " "`;
        CHROMOSOME_INDICES_SORTED+=(${elem});
done < ${SORTED_CHROMOSOMES};

# Take the nth and n - 1
CNT=${#CHROMOSOME_INDICES[@]} 
let "MAX=${CNT}/2"

declare -a methCallProcIds=()
for n in `seq 0 $MAX`; do
	let "f=${n}";
	let "r=${CNT}-1-${n}"

	if [ ! $f -gt $r ]; then
		export forward_index=${f}
		export reverse_index=${r}
	
		export CHR_INDEX_FRONT=${CHROMOSOME_INDICES_SORTED[${forward_index}]};
		export CHR_INDEX_BACK=${CHROMOSOME_INDICES_SORTED[${reverse_index}]};
		export FILE_FRONT=${FILENAME_PREFIX_METHYLATION_CALLING}.${CHR_INDEX_FRONT};
		export FILE_BACK=${FILENAME_PREFIX_METHYLATION_CALLING}.${CHR_INDEX_BACK};
		export FILE_METRICS_FRONT=${FILENAME_PREFIX_METHYLATION_CALLING_METRICS}.${CHR_INDEX_FRONT}.metrics.csv;
		export FILE_METRICS_BACK=${FILENAME_PREFIX_METHYLATION_CALLING_METRICS}.${CHR_INDEX_BACK}.metrics.csv;
		export FILE_FRONT_CP=${meth_call_files_dirname}/.${meth_call_files_basename}.${CHR_INDEX_FRONT}.checkpoint;
		export FILE_BACK_CP=${meth_call_files_dirname}/.${meth_call_files_basename}.${CHR_INDEX_BACK}.checkpoint;
	
		if [[ ! -f ${FILE_FRONT_CP} ]] || [[ ! -f ${FILE_BACK_CP} ]]; then
			if [ ! ${forward_index} -eq ${reverse_index} ]; then
				# load methylation calling script twice with a wait after another.
				FILENAME_MERGED_BAM=${bamfile} \
    				PARM_CHR_INDEX=${CHR_INDEX_FRONT} \
    				FILENAME_METH_CALLS_CHECKPOINT=${FILE_FRONT_CP} \
    				FILENAME_METH_CALLS=${FILE_FRONT} \
    				FILENAME_METH_CALL_METRICS=${FILE_METRICS_FRONT} \
    				bash ${TOOL_METHYLATION_CALLING} &
    			methCallProcIds+=($!)

				FILENAME_MERGED_BAM=${bamfile} \
    				PARM_CHR_INDEX=${CHR_INDEX_BACK} \
       				FILENAME_METH_CALLS_CHECKPOINT=${FILE_BACK_CP} \
    				FILENAME_METH_CALLS=${FILE_BACK} \
    				FILENAME_METH_CALL_METRICS=${FILE_METRICS_BACK} \
	    			bash ${TOOL_METHYLATION_CALLING} &
	    		methCallProcIds+=($!)
			else
				FILENAME_MERGED_BAM=${bamfile} \
    				PARM_CHR_INDEX=${CHR_INDEX_FRONT} \
    				FILENAME_METH_CALLS_CHECKPOINT=${FILE_FRONT_CP} \
    				FILENAME_METH_CALLS=${FILE_FRONT} \
    				FILENAME_METH_CALL_METRICS=${FILE_METRICS_FRONT} \
    				bash ${TOOL_METHYLATION_CALLING} &
    			methCallProcIds+=($!)
			fi;
		else
		    echo "Both files exist, skipping jobs for ${CHR_INDEX_FRONT}, ${CHR_INDEX_BACK}: ${FILE_FRONT}, ${FILE_BACK}";
		fi;
	fi;
done;

# wait for the methylation calling processes to be done and check if all 
# checkpoint files were created
wait ${methCallProcIds[@]} || throw $? "Error during methylation calling"

all_checkpoints_created=0
for chrom in ${CHROMOSOME_INDICES[@]}; do
	filename_meth_call_cp=${meth_call_files_dirname}/.${meth_call_files_basename}.${chrom}.checkpoint;
	if [ ! -f ${filename_meth_call_cp} ]; then
		all_checkpoints_created=1;
	fi;
done;

if [[ ${all_checkpoints_created} -eq 0 ]]; then
	touch $FILENAME_METH_CALLS_META_CHECKPOINT;
	touch $FILENAME_METH_CALLS_META_METRICS_CHECKPOINT;
fi
