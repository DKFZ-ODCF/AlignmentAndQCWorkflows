#!/bin/bash

DIR_SCRATCH=$RODDY_SCRATCH


#otherError="0"
# if the sequencing protocol is mate_pair the fastq data must be reverse complemented
# therefore "fastx_reverse_complement -Q 33" is introduced into the pipe

useReverseStepping=true
[ -z "${TRIM_END+x}" ] && TRIM_END=0 && useReverseStepping=false
[ -z "${TRIM_STEP+x}" ] && TRIM_STEP=""
[ -z "${REVERSE_STEP+x}" ] && REVERSE_STEP=""

if [ -n "${RAW_SEQ}" ]
then
    TEST_FILE=${RAW_SEQ}
    source ${TOOL_COMPRESSION_DETECTION}
    qualityScore=`${UNZIPTOOL} 2>/dev/null ${UNZIPTOOL_OPTIONS} ${RAW_SEQ}  | ${PERL_BINARY} ${TOOL_SEQUENCER_DETECTION}`

    illuminaString=""

    [ "${qualityScore}" = "phred" ] && Q=" -Q 33 "
    [ "${qualityScore}" = "illumina" ] && illuminaString="-I"

    [ "${SEQUENCER_PROTOCOL}" = "mate_pair" ] && REVERSE_STEP=" | fastx_reverse_complement ${Q} "
    # to trim the end of all readsl
    [ "${TRIM_END}" -gt 0 ] && TRIM_STEP=" | fastx_trimmer -l ${TRIM_END} ${Q} "
fi

[[ "x" == ${useAdaptorTrimming-x} ]] && useAdaptorTrimming=false # If adaptor trimming will be used this var is set. If not it is set to false

i1=""
i2=""
o1=""
o2=""
u1=""
u2=""

if [[ $useAdaptorTrimming == true ]]
then
    i1=$DIR_SCRATCH/at_i1
    i2=$DIR_SCRATCH/at_i2
    o1=$DIR_SCRATCH/at_o1
    o2=/dev/null
    u1=/dev/null
    u2=/dev/null

    mkfifo $i1 $i2 $o1
fi
