#!/bin/sh
####################################################
#
# common bwa error checking
#
####################################################

[[ -z ${TMP_FILE-} ]] && throw 11 "The variable TMP_FILE has to be set for bwa error checking"
[[ -z ${FILENAME_BWA_LOG-} ]] && throw 11 "The variable FILENAME_BWA_LOG has to be set for bwa error checking"

useMBufferStreaming="${useMBufferStreaming-false}"

# Test for success before renaming!
failure=$(grepIgnoreEmpty " fault" "$FILENAME_BWA_LOG")
[[ -n "$failure" ]] && throw 31 "found segfault '$failure' in bwa logfile!"

# if one of the read-files is terminated early (e.g. because of pipe error or error in file)
failure=$(grepIgnoreEmpty "file has fewer sequences." ${FILENAME_BWA_LOG} || true)
[[ -n "$failure" ]] && throw 41 "found segfault '$failure' in bwa logfile!"

# Barbara Aug 10 2015: I can't remember what bwa aln and sampe reported as "error".
# bluebee bwa has "error_count" in bwa-0.7.8-r2.05; and new in bwa-0.7.8-r2.06: "WARNING:top_bs_ke_be_hw: dummy be execution, only setting error."
# these are not errors that would lead to fail, in contrast to "ERROR: Bus error" 
failure=$(grepIgnoreEmpty -i "error" "$FILENAME_BWA_LOG" | grepIgnoreEmpty -v "error_count" | grepIgnoreEmpty -v "dummy be execution")
[[ -n "$failure" ]] && throw 36 "found error '$failure' in bwa logfile!"

failure=$(grepIgnoreEmpty "Abort. Sorry." "$FILENAME_BWA_LOG")
[[ -n "$failure" ]] && throw 37 "found error '$failure' in bwa logfile!"

# An assumption of the workflow is violated.
failure=$(grepIgnoreEmpty "file has fewer sequences." "$FILENAME_BWA_LOG")
[[ -n "$failure" ]] && throw 41 "found error '$failure' in bwa logfile!"

# samtools sort may complain about truncated temp files and for each line outputs
# the error message. This happens when the same files are written at the same time,
# see http://sourceforge.net/p/samtools/mailman/samtools-help/thread/BAA90EF6FE3B4D45A7B2F6E0EC5A8366DA3AB5@USTLMLLYC102.rf.lilly.com/
# This happens when the scheduler puts the same job on 2 nodes bc. the prefix for samtools-0.1.19 -o $prefix is constructed using the job ID
echo "testing $FILENAME_SORT_LOG"
if [[ -n "$FILENAME_SORT_LOG" && -f $FILENAME_SORT_LOG ]]; then
	failure=$(grepIgnoreEmpty "is truncated. Continue anyway." "$FILENAME_SORT_LOG")
	if [[ -n "$failure" ]]; then
	    throw 38 "found error '$failure' in samtools sorting logfile!"
	fi
else
	echo "there is no samtools sort log file"
fi

if [[ $(cat "$FILENAME_BWA_EC" | cut -f 1 -d ' ') != "0" ]]; then
    throw 32 "There was a non-zero exit code in bwa pipe (w/ pipefail): $(cat \"$FILENAME_BWA_EC\"); exiting..."
fi

# disable file size checks if the interprocess streaming is active
if [[ $useMBufferStreaming != true ]]; then

    if [[ "2048" -gt `stat -c %s "$TMP_FILE"` ]]; then
        errorString="Output file is too small!"
        throw 33 "$errorString"
    fi

    # TODO Check that?
    if [[ "$?" != "0"  ]]; then
        throw 32 "$errorString"
    fi
fi

echo "all OK"
