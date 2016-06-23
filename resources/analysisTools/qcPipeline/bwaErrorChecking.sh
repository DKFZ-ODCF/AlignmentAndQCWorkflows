#!/bin/sh
####################################################
#
# common bwa error checking
#
####################################################

[[ -z ${TMP_FILE-} ]] && echo "The variable TMP_FILE has to be set for bwa error checking" && exit 11
[[ -z ${FILENAME_BWA_LOG-} ]] && echo "The variable FILENAME_BWA_LOG has to be set for bwa error checking" && exit 11

useMBufferStreaming=${useMBufferStreaming-false}

# disable file size checks if the interprocess streaming is active
if [[ $useMBufferStreaming != true ]]
then

    if [[ "2048" -gt `stat -c %s ${TMP_FILE}` ]]
    then
        errorString="Output file is too small!"
        exit 33
    fi

    # TODO Check that?
    [[ "$?" != "0"  ]] && echo $errorString && exit 32
fi

# test for success before renaming!
success=`grep " fault" ${FILENAME_BWA_LOG}`
[[ ! -z "$success" ]] && echo found segfault $success in bwa logfile! && exit 31

success=`grep -i "error" ${FILENAME_BWA_LOG} | grep -v "error_count"`
[[ ! -z "$success" ]] && echo found error $success in bwa logfile! && exit 36


# samtools sort may complain about truncated temp files and for each line outputs
# the error message. This happens when the same files are written at the same time,
# see http://sourceforge.net/p/samtools/mailman/samtools-help/thread/BAA90EF6FE3B4D45A7B2F6E0EC5A8366DA3AB5@USTLMLLYC102.rf.lilly.com/
# This happens when the scheduler puts the same job on 2 nodes bc. the prefix for samtools-0.1.19 -o $prefix is constructed using the job ID
echo "testing $FILENAME_SORT_LOG"
if [ ! -z $FILENAME_SORT_LOG ] && [ -f $FILENAME_SORT_LOG ]
then
	success=`grep "is truncated. Continue anyway." $FILENAME_SORT_LOG`
	[[ ! -z "$success" ]] && echo found error $success in samtools sorting logfile! && exit 36
else
	echo "there is no samtools sort log file"
fi
echo all OK
