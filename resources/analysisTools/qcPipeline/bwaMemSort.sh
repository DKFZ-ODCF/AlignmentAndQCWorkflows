#!/bin/bash

source "$TOOL_BASH_LIB"

printInfo

set -o pipefail

ON_CONVEY=${useAcceleratedHardware:-false}

# use scratch dir for temp files: samtools sort uses the current working directory for them
WORKDIR=${DIR_TEMP}/${RODDY_JOBID}
LCLSCRATCH=${RODDY_SCRATCH}
mkdir -p ${WORKDIR}
cd ${WORKDIR}

FNPIPE1=${LCLSCRATCH}/NAMED_PIPE1
FNPIPE2=${LCLSCRATCH}/NAMED_PIPE2
FNPIPE_FLAGSTATS_IN=${LCLSCRATCH}/NP_FLAGSTATS_IN
FNPIPE_INDEX_IN=${LCLSCRATCH}/NP_INDEX_IN

#NP_FLAGSTATS=$RODDY_SCRATCH/NAMED_PIPE_FLAGSTATS
#NP_INDEX=$RODDY_SCRATCH/NAMED_PIPE_FLAGSTATS
#NP_ISIZES=$RODDY_SCRATCH/NAMED_PIPE_ISIZES

MBUFFER_2="mbuffer -q -m ${MEMSORT_MBUFFER_SIZE-8}G -l /dev/null"

# the more efficient version pipes via local scratch dir (to avoid network problems)
mkfifo ${FNPIPE1}
mkfifo ${FNPIPE2}
mkfifo ${FNPIPE_FLAGSTATS_IN}
mkfifo ${FNPIPE_INDEX_IN}
echo created pipes

INDEX_FILE=${FILENAME_SORTED_BAM}.bai

#Samtools always attaches .bam, so we use a different output filename for samtools without .bam
#TMP_FILE_INDEX=${FILENAME_SORTED_BAM}.bai
TMP_FILE_SAMTOOLS=${FILENAME_SORTED_BAM}_tmp
TMP_FILE=${TMP_FILE_SAMTOOLS}
TMP_FILE_SORT=${WORKDIR}/sorted.temporary
TMP_FILE_INDEX=${FILENAME_SORTED_BAM}_tmp.bai
FLAGSTAT_TMP=${FILENAME_FLAGSTAT}_tmp


# error tracking!
FILENAME_BWA_LOG=${DIR_TEMP}/`basename ${FILENAME_SORTED_BAM}`_errlog_sampesort
#LOG_INDEX=${FILENAME_SORTED_BAM}_errlog_index

RAW_SEQ=${RAW_SEQ_1}
source ${TOOL_COMMON_ALIGNMENT_SETTINGS_SCRIPT}


#source ${TOOLSDIR}/../tools/createProcessUtilizationratesLogfile.sh

SAMTOOLS_SORT_BINARY=samtools-0.1.19
#
#if [[ ! -s ${RAW_SEQ_2} ]] # file size of second file is zero
#then
#    # Single end read is active
#    nice ${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_1 ${TRIM_STEP} ${REVERSE_STEP} | mbuffer -q -m 1G -l /dev/null > $FNPIPE1 &
#    ${BWA_BINARY} samse -r "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" ${INDEX_PREFIX} ${FILENAME_SEQ_1} $FNPIPE1 2> ${FILENAME_BWA_LOG} | $MBUFFER_2 | ${SAMTOOLS_BINARY} view -uSbh - | $MBUFFER_2 | ${SAMTOOLS_SORT_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} -o - ${TMP_FILE_SORT} | tee ${TMP_FILE_SAMTOOLS} | ${SAMTOOLS_SORT_BINARY} flagstat - > ${FLAGSTAT_TMP}

#if [[ $useAdaptorTrimming == false ]]
#then
#
#    ##Reset error and tmp variables. these are modified by bwaCommonAlignment...
#    nice ${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_1 ${TRIM_STEP} ${REVERSE_STEP} | $MBUFFER_2 > $FNPIPE1 &
#
#    # Repeat some steps for the second named pipe
#    nice ${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} $RAW_SEQ_2 ${TRIM_STEP} ${REVERSE_STEP} | $MBUFFER_2 > $FNPIPE2 &
#
#
#    ## -r STR   read group header line such as `@RG\tID:foo\tSM:bar' [null]
#    #Samtools 0.1.19 uses multithreaded sorting => therefore it is hardcoded here.
#    ${BWA_BINARY} mem -t ${THREADS} -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" $BWA_MEM_OPTIONS ${INDEX_PREFIX} $FNPIPE1 $FNPIPE2 | $MBUFFER_2 | ${SAMTOOLS_BINARY} view -uSbh - | $MBUFFER_2 | ${SAMTOOLS_SORT_BINARY} sort -@ ${THREADS} -m ${SAMPESORT_MEMSIZE} -o - ${TMP_FILE_SORT} | tee ${TMP_FILE_SAMTOOLS} | ${SAMTOOLS_SORT_BINARY} flagstat - > ${FLAGSTAT_TMP}
##    ${BWA_BINARY} sampe -P -T -t 8 ${BWA_SAMPESORT_OPTIONS} -r "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" ${INDEX_PREFIX} ${FILENAME_SEQ_1} ${FILENAME_SEQ_2} $FNPIPE1 $FNPIPE2 2> ${FILENAME_BWA_LOG} | $MBUFFER_2 | ${SAMTOOLS_BINARY} view -uSbh - | $MBUFFER_2 | ${SAMTOOLS_SORT_BINARY} sort -@ 8 -m ${SAMPESORT_MEMSIZE} -o - ${TMP_FILE_SORT} | tee ${TMP_FILE_SAMTOOLS} | ${SAMTOOLS_SORT_BINARY} flagstat - > ${FLAGSTAT_TMP}
#
#else

if [[ $useAdaptorTrimming == true ]] # [[ "$ADAPTOR_TRIMMING_TOOL" == *.jar ]]
then
    eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | $MBUFFER_2 > $i1" &
    eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | $MBUFFER_2 > $i2" &

    #We need a second output pipe for sampe
    o2=${WORKDIR}/at_o2
    mkfifo $o2
    eval "$JAVA_BINARY -jar ${TOOL_ADAPTOR_TRIMMING} $ADAPTOR_TRIMMING_OPTIONS_0 $i1 $i2 $o1 $u1 $o2 $u2 $ADAPTOR_TRIMMING_OPTIONS_1" &
    cat $o1 ${TRIM_STEP} ${REVERSE_STEP} | $MBUFFER_2 > $FNPIPE1 &
    cat $o2 ${TRIM_STEP} ${REVERSE_STEP} | $MBUFFER_2 > $FNPIPE2 &
else
    eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | $MBUFFER_2 > $FNPIPE1" &
    eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | $MBUFFER_2 > $FNPIPE2" &
fi
useBioBamBamSort=${useBioBamBamSort-true}

NP_BAMSORT=${LCLSCRATCH}/NAMED_PIPE_BAMSORT
FILENAME_BAMSORT_TMP=${WORKDIR}/`basename ${FILENAME_SORTED_BAM}`.tmp
mkfifo ${NP_BAMSORT}
LENGTH_SEQ_1=`${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | head | wc -l`
LENGTH_SEQ_2=`${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | head | wc -l`

INPUT_PIPES=""
[[ ${LENGTH_SEQ_1} != 0 ]] && INPUT_PIPES="$FNPIPE1"
[[ ${LENGTH_SEQ_2} != 0 ]] && INPUT_PIPES="${INPUT_PIPES} $FNPIPE2"
[[ ${LENGTH_SEQ_1} == 0 ]] && cat $FNPIPE1 >/dev/null
[[ ${LENGTH_SEQ_2} == 0 ]] && cat $FNPIPE2 >/dev/null


if [[ "$ON_CONVEY" == "true" ]]; then
    bamname=`basename ${FILENAME_SORTED_BAM}`
    useBioBamBamSort=false;
    bwaBinary=${BWA_ACCELERATED_BINARY}
    cat ${FNPIPE_FLAGSTATS_IN} | ${SAMBAMBA_BINARY} flagstat /dev/stdin > $FLAGSTAT_TMP & procID_FSTATS=$!
    (set -o pipefail; ${bwaBinary} mem -t 12 -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" $BWA_MEM_OPTIONS ${INDEX_PREFIX} ${INPUT_PIPES} | $MBUFFER_2 | \
        ${SAMBAMBA_BINARY} view -f bam -S -l 0 -t 8 /dev/stdin | $MBUFFER_2 | \
        tee ${FNPIPE_FLAGSTATS_IN} | ${SAMBAMBA_BINARY} sort --tmpdir=${LCLSCRATCH} -l 9 -t 8 -m ${SAMPESORT_MEMSIZE} /dev/stdin -o $TMP_FILE; echo $? > ${DIR_TEMP}/${bamname}_ec) & procID_MEMSORT=$!
    wait $procID_MEMSORT
    [[ ! `cat ${DIR_TEMP}/${bamname}_ec` -eq "0" ]] && echo "bwa mem-sambamba pipe returned a non-zero exit code and the job will die now." && exit 100
    wait $procID_FSTATS
    [[ ! $? -eq 0 ]] && exit 100
elif [[ ${useBioBamBamSort} == false ]]; then
    cat ${FNPIPE_INDEX_IN} | ${SAMTOOLS_SORT_BINARY} index - ${TMP_FILE_INDEX} & procID_IDX=$!
    (set -o pipefail; ${BWA_BINARY} mem -t ${THREADS} -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" $BWA_MEM_OPTIONS ${INDEX_PREFIX} ${INPUT_PIPES} | $MBUFFER_2 | \
        ${SAMTOOLS_BINARY} view -uSbh - | $MBUFFER_2 | \
        ${SAMTOOLS_SORT_BINARY} sort -@ ${THREADS} -m ${SAMPESORT_MEMSIZE} -o - ${TMP_FILE_SORT} | tee ${TMP_FILE_SAMTOOLS} ${FNPIPE_INDEX_IN} | ${SAMTOOLS_SORT_BINARY} flagstat - > ${FLAGSTAT_TMP}; echo $? > ${DIR_TEMP}/${bamname}_ec) & procID_MEMSORT=$!
    wait $procID_MEMSORT
    [[ ! `cat ${DIR_TEMP}/${bamname}_ec` -eq "0" ]] && echo "bwa mem-samtools pipe returned a non-zero exit code and the job will die now." && exit 100
    wait procID_IDX
    [[ ! $? -eq 0 ]] && exit 100
else
    (set -o pipefail; ${BWA_BINARY} mem -t ${THREADS} -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" $BWA_MEM_OPTIONS ${INDEX_PREFIX} ${INPUT_PIPES} | $MBUFFER_2 | \
        ${SAMTOOLS_BINARY} view -uSbh - | $MBUFFER_2 | \
        ${BAMSORT_BINARY} \
            O=${NP_BAMSORT} \
            level=1 \
            inputthreads=2 \
            outputthreads=2 \
            index=1 \
            indexfilename=${INDEX_FILE} \
            calmdnm=1 \
            calmdnmrecompindetonly=1 \
            calmdnmreference=${INDEX_PREFIX} \
            tmpfile=${FILENAME_BAMSORT_TMP}; echo $? > ${DIR_TEMP}/ec_bbam) &

    cat ${NP_BAMSORT} | tee ${TMP_FILE_SAMTOOLS} | ${SAMTOOLS_SORT_BINARY} flagstat - > ${FLAGSTAT_TMP}
    wait
    sleep 30
fi

rm $i1 $i2 $o1 $o2 2> /dev/null
rm -rf ${FILENAME_BAMSORT_TMP}* 2> /dev/null # Remove temporary files.

#fi
#${SAMTOOLS_BINARY} index ${NP_INDEX}
errorString="There was a non-zero exit code in the bwa sampe - samtools sort pipeline; exiting..." 

[[ -f ${TMP_FILE_INDEX} ]] && mv ${TMP_FILE_INDEX} ${INDEX_FILE}

source ${TOOL_BWA_ERROR_CHECKING_SCRIPT}
[[ "$ON_CONVEY" == "false" ]] && [[ ! -s ${INDEX_FILE} ]] && echo "Bam Index is of size 0; Exitting" && exit 5

#rm $FILE_BENCHMARK_STAYALIVE
#wait $utilizationLogProcess

mv ${FLAGSTAT_TMP} ${FILENAME_FLAGSTAT}
mv ${TMP_FILE_SAMTOOLS} ${FILENAME_SORTED_BAM}
[[ "$ON_CONVEY" == "false" ]] && touch ${INDEX_FILE}


rm $FNPIPE1
rm $FNPIPE2
