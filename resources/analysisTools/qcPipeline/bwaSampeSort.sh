#!/bin/bash

#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=6
#PBS -l mem=52g
#PBS -m a

source ${CONFIG_FILE}

set -o pipefail

# use scratch dir for temp files: samtools sort uses the current working directory for them
LCLSCRATCH=${PBS_SCRATCH_DIR}/${PBS_JOBID}
cd $LCLSCRATCH

FNPIPE1=$LCLSCRATCH/NAMED_PIPE1
FNPIPE2=$LCLSCRATCH/NAMED_PIPE2
FNPIPE_FLAGSTATS=${LCLSCRATCH}/FNPIPE_FLAGSTATS
FNPIPE_INDEX=${LCLSCRATCH}/FNPIPE_INDEX

MBUFFER_2="mbuffer -q -m 8G -l /dev/null"

# the more efficient version pipes via local scratch dir (to avoid network problems)
mkfifo $FNPIPE1 $FNPIPE2 $FNPIPE_FLAGSTATS $FNPIPE_INDEX
echo created pipes

#Samtools always attaches .bam, so we use a different output filename for samtools without .bam
#TMP_FILE_INDEX=${FILENAME_SORTED_BAM}.bai
TMP_FILE_SAMTOOLS=${FILENAME_SORTED_BAM}_tmp
TMP_FILE=${TMP_FILE_SAMTOOLS}       # This variable is used for the error checking code.
FILENAME_INDEX=${FILENAME_SORTED_BAM}.bai
TMP_FILE_SORT=${LCLSCRATCH}/sorted.temporary
TMP_FILE_FLAGSTATS=${FILENAME_FLAGSTAT}_tmp
TMP_FILE_INDEX=${FILENAME_INDEX}.tmp

# error tracking!
FILENAME_BWA_LOG=${DIR_TEMP}/`basename ${FILENAME_SORTED_BAM}`_errlog_sampesort

RAW_SEQ=${RAW_SEQ_1}
source ${TOOL_COMMON_ALIGNMENT_SETTINGS_SCRIPT}

SAMTOOLS_SORT_BINARY=samtools-0.1.19

bwaBinary="${BWA_BINARY}"
useConvey=false
# If things go wrong with this, one could use which with the convey binary!
if [[ ${PBS_QUEUE} == "convey" ]]
then
    bwaBinary=${BWA_ACCELERATED_BINARY}
    useConvey=true
fi

${SAMTOOLS_BINARY} index $FNPIPE_INDEX $TMP_FILE_INDEX &
cat $FNPIPE_FLAGSTATS | ${SAMTOOLS_SORT_BINARY} flagstat - > ${TMP_FILE_FLAGSTATS} &

useSingleEndProcessing=${useSingleEndProcessing-false}

if [[ ${useSingleEndProcessing} == "true" ]]
then
    # Single end read is active
    nice ${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_1 ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 1G -l /dev/null > $FNPIPE1 &
    ${BWA_BINARY} samse -r "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" ${INDEX_PREFIX} ${FILENAME_SEQ_1} $FNPIPE1 2> ${FILENAME_BWA_LOG} | $MBUFFER_2 | ${SAMTOOLS_BINARY} view -uSbh - | $MBUFFER_2 | ${SAMTOOLS_SORT_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} -o - ${TMP_FILE_SORT} | tee ${TMP_FILE_SAMTOOLS} $FNPIPE_INDEX $FNPIPE_FLAGSTATS > /dev/null
else
    if [[ $useAdaptorTrimming == false ]]
    then

        ##Reset error and tmp variables. these are modified by bwaCommonAlignment...
        eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_1 ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 1G -l /dev/null > ${FNPIPE1}" &

        # Repeat some steps for the second named pipe
        eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_2 ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 1G -l /dev/null > ${FNPIPE2}" &

    else
        eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | mbuffer -q -m 6G -l /dev/null > $i1" &
        eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | mbuffer -q -m 6G -l /dev/null > $i2" &

        #We need a second output pipe for sampe
        o2=${DIR_SCRATCH}/at_o2
        mkfifo $o2

        if [[ "$ADAPTOR_TRIMMING_TOOL" == *.jar ]]
        then
            eval "java7 -jar  ${TOOL_ADAPTOR_TRIMMING} $ADAPTOR_TRIMMING_OPTIONS_0 $i1 $i2 $o1 $u1 $o2 $u2 $ADAPTOR_TRIMMING_OPTIONS_1" &
        fi

        cat $o1 ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 6G -l /dev/null > $FNPIPE1 &
        cat $o2 ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 6G -l /dev/null > $FNPIPE2 &

    fi
    ${BWA_BINARY} sampe -P -T -t 8 ${BWA_SAMPESORT_OPTIONS} -r "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" ${INDEX_PREFIX} ${FILENAME_SEQ_1} ${FILENAME_SEQ_2} $FNPIPE1 $FNPIPE2 2> ${FILENAME_BWA_LOG} | $MBUFFER_2 | ${SAMTOOLS_BINARY} view -uSbh - | $MBUFFER_2 | ${SAMTOOLS_SORT_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} -o - ${TMP_FILE_SORT} | tee ${TMP_FILE_SAMTOOLS} $FNPIPE_INDEX $FNPIPE_FLAGSTATS > /dev/null
fi

#Wait for unfinished samtools jobs.
wait

#${SAMTOOLS_BINARY} index ${NP_INDEX}
errorString="There was a non-zero exit code in the bwa sampe - samtools sort pipeline; exiting..." 

source ${TOOL_BWA_ERROR_CHECKING_SCRIPT}

#rm $FILE_BENCHMARK_STAYALIVE
#wait $utilizationLogProcess

mv ${TMP_FILE_FLAGSTATS} ${FILENAME_FLAGSTAT}
mv ${TMP_FILE_INDEX} ${FILENAME_INDEX}
mv ${TMP_FILE_SAMTOOLS} ${FILENAME_SORTED_BAM}

rm $FNPIPE1
rm $FNPIPE2
