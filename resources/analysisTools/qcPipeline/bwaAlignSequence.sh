#!/bin/bash

#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=12
#PBS -S /bin/bash
#PBS -l mem=3600m
#PBS -m a

source "$TOOL_BASH_LIB"

printInfo

set -o pipefail

ON_CONVEY=$(runningOnConvey)

TMP_FILE=${FILENAME_ALIGNMENT}_temp
# error tracking because bwa never exits anything but 0
FILENAME_BWA_LOG=${DIR_TEMP}/`basename ${FILENAME_ALIGNMENT}`_errlog

source ${TOOL_COMMON_ALIGNMENT_SETTINGS_SCRIPT}

# Check for single lane processing.
# For single lanes you also need two raw sequences but one must be of size 0.
# If the lane for this job is of size 0 the aligned file is created with size 0.
# Afterwards the job exits with exit code 0
[[ -z "${useSingleEndProcessing+x}" ]] && useSingleEndProcessing=false

if [[ "${useSingleEndProcessing}" == "true" && ! -f ${RAW_SEQ} ]]
then
    touch ${FILENAME_ALIGNMENT}
    exit 0
fi

cmd=""
bwaBinary="${BWA_BINARY}"
alnThreadOptions="8"
useConvey=false

# If things go wrong with this, one could use which with the convey binary!
if [[ "$ON_CONVEY" == "true" ]]
then
    bwaBinary=${BWA_ACCELERATED_BINARY}
    useConvey=true
    alnThreadOptions="12"
fi

# Check the bwa version
sh ${TOOL_CHECK_BWA_AND_INDEX_VERSIONS} ${bwaBinary} ${INDEX_PREFIX}
[[ $? -ne 0 ]] && echo "Problems when checking bwa binary against the index prefix." && exit -5
# [[ ${useReverseStepping} ]] &&

baseBWACall="${bwaBinary} aln -t ${alnThreadOptions} ${BWA_ALIGNMENT_OPTIONS} ${illuminaString} ${INDEX_PREFIX}"



#if [[ ${useConvey} = true && ${useReverseStepping} = false ]]
#then
#    cmd="${baseBWACall} ${PRM_RAW_SEQ} > ${TMP_FILE} 2> ${FILENAME_BWA_LOG}"
#else
    # In this case either the convey with reverse stepping is used or the default queue with or without reverse stepping.
#    cmd="${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${PRM_RAW_SEQ} ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 100M -l /dev/null | ${baseBWACall} - > ${TMP_FILE} 2> ${FILENAME_BWA_LOG}"
#fi

#source ${TOOLSDIR}/../tools/createProcessUtilizationratesLogfile.sh

if [[ $useAdaptorTrimming == false ]]
then
    eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ} ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 100M -l /dev/null | ${baseBWACall} - > ${TMP_FILE} 2> ${FILENAME_BWA_LOG}"
else
    eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ}   | mbuffer -q -m 2G  -l /dev/null > $i1" &
    eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | mbuffer -q -m 2G  -l /dev/null > $i2" &

    if [[ "$RAW_SEQ_FILE_1_INDEX" -lt "$RAW_SEQ_FILE_2_INDEX" ]]
    then
        if [[ ${TOOL_ADAPTOR_TRIMMING} == *.jar ]]
        then
            eval "java7 -jar  ${TOOL_ADAPTOR_TRIMMING} $ADAPTOR_TRIMMING_OPTIONS_0 $i1 $i2 $o1 $u1 $o2 $u2 $ADAPTOR_TRIMMING_OPTIONS_1" &
        fi

    else
        # Roddy calls alignment on a per lane file base. That means that the process is called the first time with PRM_RAW_SEQ having index = 1 and the second time with PRM_RAW_SEQ having index = 2
        # In the second case, the input files for the adaptor trimming with trimmomatic need to be swapped, as trimmomatic needs the forward read as file1 and the reverse read as file2.
        # trimmomatic is therefore called with i2 i1 as the input. Also the output o1 and o2 are swapped to match the input.
        if [[ ${TOOL_ADAPTOR_TRIMMING} == *.jar ]]
        then
            eval "java7 -jar  ${TOOL_ADAPTOR_TRIMMING} $ADAPTOR_TRIMMING_OPTIONS_0 $i2 $i1 $o2 $u1 $o1 $u2 $ADAPTOR_TRIMMING_OPTIONS_1" &
        fi

    fi

    eval "cat $o1 ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 2G -l /dev/null | ${baseBWACall} - | mbuffer -m 2G > ${TMP_FILE} 2> ${FILENAME_BWA_LOG}"

    rm $i1 $i2 $o1
fi


errorString="There was a non-zero exit code in bwa aln; exiting..." 
source ${TOOL_BWA_ERROR_CHECKING_SCRIPT}

#rm $FILE_BENCHMARK_STAYALIVE
#wait $utilizationLogProcess

mv ${TMP_FILE} ${FILENAME_ALIGNMENT}

