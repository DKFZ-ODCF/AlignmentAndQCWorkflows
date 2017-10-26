#!/bin/bash
#PBS -l walltime=0:10:00
#PBS -l nodes=1
#PBS -m a

# Pipe in QC_summary.txt from analysis dir...

##job writeQCsummary_newFlagstat
cmd=" -p ${PID} -s ${SAMPLE} -r ${RUN} -l ${LANE} -w ${FILENAME_QCSUM}_WARNING.txt"	# new!
[[ -n ${FILENAME_FLAGSTAT+x} ]] && cmd=$cmd" -f ${FILENAME_FLAGSTAT}"
[[ -n ${FILENAME_DIFFCHROM+x} ]] && cmd=$cmd" -d ${FILENAME_DIFFCHROM}"
[[ -n ${FILENAME_ISIZE+x} ]] && cmd=$cmd" -i ${FILENAME_ISIZE}"
[[ -n ${FILENAME_COVERAGE+x} ]] && cmd=$cmd" -c ${FILENAME_COVERAGE}"
[[ -n ${FILENAME_TARGETCAPTURE+x} ]] && cmd=$cmd" -t ${FILENAME_TARGETCAPTURE}"
[[ -n ${FILENAME_METRICS+x} ]] && cmd=$cmd" -m ${FILENAME_METRICS}"


${PERL_BINARY} ${TOOL_WRITE_QC_SUMMARY} $cmd > ${FILENAME_QCSUM}
