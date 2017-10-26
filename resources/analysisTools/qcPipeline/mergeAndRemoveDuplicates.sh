#!/bin/bash

#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=6
#PBS -l mem=50g
#PBS -m a

set -o pipefail

today=`date +'%Y-%m-%d_%Hh%M'`
WORKDIR=${DIR_TEMP}/${RODDY_JOBID}
LCLSCRATCH=${RODDY_SCRATCH}
TMP_DIR=${FILENAME}_PICARDTEMP/
mkdir $TMP_DIR
mkdir -p $WORKDIR

SAM_BIN_MT=samtools-0.1.19

TMP_FILE=${TMP_DIR}FILENAME_TMP
TMP_FLAGSTATS=${TMP_FILE}_flagstats
TMP_INDEX=${TMP_FILE}.bai.tmp
NP_PIC_OUT=${LCLSCRATCH}/np_picard_out.sam
NP_COMP_IN=${LCLSCRATCH}/np_compression_in
NP_METRICS_OUT=${LCLSCRATCH}/np_metrics_out
NP_INDEX_IN=${LCLSCRATCH}/np_index_in
NP_FLAGSTATS_IN=${LCLSCRATCH}/np_flagstats_in
filenameRCMarkDup=${TMP_DIR}/rcMarkDup.txt

# if the merged file already exists, only merge new lanes to it
if [ -e ${FILENAME} ] && [ -s ${FILENAME} ]
then
        singlebams=`${SAMTOOLS_BINARY} view -H ${FILENAME} | grep "^@RG"`
        [[ -z "$singlebams" ]] && echo "could not detect single lane BAM files in ${FILENAME}, stopping here" && exit 23

        notyetmerged=`perl ${TOOL_CHECK_ALREADY_MERGED_LANES} ${INPUT_FILES} "$singlebams" ${pairedBamSuffix}`
        [[ "$?" != 0 ]] && echo "something went wrong with the detection of merged files in ${FILENAME}, stopping here" && exit 24

        # the Perl script returns BAM names separated by :, ending with :
        if [ -z $notyetmerged ]
        then
                echo "All listed BAM files are already in ${FILENAME}, not running merge-markdup again."
                echo "Checking flagstats and index files for existence, if inexistent they are created again."
                ([[ ! -f ${FILENAME_FLAGSTATS} ]] && ${SAM_BIN_MT} flagstat ${FILENAME} > ${TMP_FLAGSTATS} && mv ${TMP_FLAGSTATS} ${FILENAME_FLAGSTATS}) & pIDFstats=$!
                ([[ ! -f ${TMP_INDEX} ]] && ${SAM_BIN_MT} index ${FILENAME} ${TMP_INDEX} && mv ${TMP_INDEX} ${FILENAME}.bai) & pIDIdx=$!

                wait $pIDFstats; [[ $? -gt 0 ]] && echo "Error in flagstats" && exit 5
                wait $pIDIdx; [[ $? -gt 0 ]] && echo "Error in index" && exit 5

                exit 0
        else
                # input files is now the merged file and the new file(s)
                INPUT_FILES=${FILENAME}":"$notyetmerged
                # keep the old metrics file for comparison
                mv ${FILENAME_METRICS} ${FILENAME_METRICS}_before_${today}.txt
        fi
fi


mkfifo ${NP_PIC_OUT} ${NP_COMP_IN} ${NP_METRICS_OUT} ${NP_INDEX_IN} ${NP_FLAGSTATS_IN}

oIFS=${IFS}
IFS=":"
mergeINPUT=''
for bamfile in ${INPUT_FILES}; do mergeINPUT=${mergeINPUT}' I='$bamfile; done
IFS=${oIFS}


MBUF_1G="mbuffer -m 1g -q -l /dev/null"
MBUF_2G="mbuffer -m 2g -q -l /dev/null"
MBUF_4G="mbuffer -m 2g -q -l /dev/null"

useBioBamBamMarkDuplicates=${useBioBamBamMarkDuplicates-true}

if [[ ${useBioBamBamMarkDuplicates} == false ]]; then
    FILEHANDLES=$((`ulimit -n` - 16))

    # Complicate error catching because:
    # In some cases picard exits before cat starts to work on the named pipes. If so, the cat processes block and wait for input data forever.
    # TODO Think hard if it helps to just move picard two steps back and to start the cat processes earlier.

    (JAVA_OPTIONS="-Xms64G -Xmx64G" picard.sh MarkDuplicates ${mergeINPUT} OUTPUT=${NP_PIC_OUT} METRICS_FILE=${NP_METRICS_OUT} MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${FILEHANDLES} TMP_DIR=${TMP_DIR} COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT ${mergeAndRemoveDuplicates_argumentList} ASSUME_SORTED=TRUE CREATE_INDEX=FALSE MAX_RECORDS_IN_RAM=12500000; echo $? > ${filenameRCMarkDup}) & procIDPicard=$!
    (cat ${NP_METRICS_OUT} | ${MBUF_1G} > ${FILENAME_METRICS}) & procIDMetricsPipe=$!
    (cat ${NP_PIC_OUT} | ${MBUF_2G}  > ${NP_COMP_IN}) & procIDPicardOutPipe=$!
    (${SAM_BIN_MT} view -S -@ 8 -b ${NP_COMP_IN} | ${MBUF_2G} | tee ${NP_INDEX_IN} ${NP_FLAGSTATS_IN} | ${MBUF_4G} > ${TMP_FILE}) & procIDSamtoolsView=$!

    (${SAM_BIN_MT} index ${NP_INDEX_IN}) & procIDSamtoolsIndex=$!
    (${SAM_BIN_MT} flagstat ${NP_FLAGSTATS_IN} > ${TMP_FLAGSTATS}) & procIDSamtoolsFlagstat=$!

    wait $procIDPicard
    [[ ! `cat ${filenameRCMarkDup}` -eq "0" ]] && echo "Picard returned an exit code and the job will die now." && exit 100

    wait $procIDMetricsPipe $procIDPicardOutPipe $procIDSamtoolsFlagstat $procIDSamtoolsIndex $procIDSamtoolsView

else
#        rewritebam=1 \
#        rewritebamlevel=1 \
    (${BAMMARKDUPLICATES_BINARY} \
        M=${FILENAME_METRICS} \
        tmpfile=${WORKDIR}/biobambammerge.tmp \
        markthreads=8 \
        level=9 \
        index=1 \
        indexfilename=${FILENAME}.bai \
        ${mergeINPUT} \
        O=${NP_PIC_OUT}; echo $? > ${filenameRCMarkDup}) &

    cat ${NP_PIC_OUT} | ${MBUF_4G} | tee ${NP_FLAGSTATS_IN} | ${MBUF_4G} > ${TMP_FILE} &
    ${SAM_BIN_MT} flagstat ${NP_FLAGSTATS_IN} > ${TMP_FLAGSTATS} &

    wait
fi


[[ "$?" != 0 ]] && echo "Non-zero exit code from picard MergeSamFiles; exiting..." && exit 2

mv ${TMP_FILE} ${FILENAME}
[[ ${useBioBamBamMarkDuplicates} == false ]] && mv ${NP_INDEX_IN}.bai ${FILENAME}.bai
[[ ${useBioBamBamMarkDuplicates} == true  ]] && [[ ! -f ${FILENAME}.bai ]] && echo "Bam Index file does not exist but should do so!" && exit 5
[[ ${useBioBamBamMarkDuplicates} == true  ]] && touch ${FILENAME}.bai   # Update timestamp

mv ${TMP_FLAGSTATS} ${FILENAME_FLAGSTATS}

rm -rf $TMP_DIR
