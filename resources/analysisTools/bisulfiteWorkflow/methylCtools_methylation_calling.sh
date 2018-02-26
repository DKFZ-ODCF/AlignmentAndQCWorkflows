#!/bin/bash

source ${TOOL_BASH_LIB}

set -o pipefail

set -xuv

# Get Parameters
CHR_PREFIX=${CHR_PREFIX-}
CHR_SUFFIX=${CHR_SUFFIX-}
FILENAME_MERGED_BAM=${FILENAME_MERGED_BAM-}
FILENAME_METH_CALLS=${FILENAME_METH_CALLS-}
FILENAME_METH_CALL_METRICS=${FILENAME_METH_CALL_METRICS-}
PARM_CHR_INDEX=${PARM_CHR_INDEX-}
FILENAME_METH_CALLS_BASE=`basename ${FILENAME_METH_CALLS}`

# Create named pipes
NP_METHCALLS_CG=${RODDY_SCRATCH}/${FILENAME_METH_CALLS_BASE}.CG
NP_METHCALLS_CH=${RODDY_SCRATCH}/${FILENAME_METH_CALLS_BASE}.CH
mkfifo ${NP_METHCALLS_CG} ${NP_METHCALLS_CH}

# Check if all parameters are set
[[ -z ${FILENAME_MERGED_BAM} ]] && throw -10 "Parameter is missing: FILENAME_MERGED_BAM"
[[ -z ${FILENAME_METH_CALLS} ]] && throw -10 "Parameter is missing: FILENAME_METH_CALLS"
[[ -z ${FILENAME_METH_CALL_METRICS} ]] && throw -10 "Parameter is missing: FILENAME_METH_CALL_METRICS"

# Construct current sequence name
chr=""
if [[ -n "$PARM_CHR_INDEX" ]]; then
    chr=${CHR_PREFIX}${PARM_CHR_INDEX}${CHR_SUFFIX}
fi

# Get directories and create filenames
METH_CALL_SUBDIR=`dirname ${FILENAME_METH_CALLS}`
SCRATCH_DIR=${RODDY_SCRATCH}

METH_CALL_METRICS_SUBDIR=`dirname ${FILENAME_METH_CALL_METRICS}`

FILENAME_METH_CALL_CG_TMP=${FILENAME_METH_CALLS}.CG.bed.tmp
FILENAME_METH_CALL_CH_TMP=${FILENAME_METH_CALLS}.CH.bed.tmp
FILENAME_METH_CALL_CG=${FILENAME_METH_CALLS}.CG.bed
FILENAME_METH_CALL_CH=${FILENAME_METH_CALLS}.CH.bed

################################################
# Conversion of Methylation Calls -- CG and CH #
################################################
if [[ ${METH_CALLS_CONVERTER} == 'moabs' ]]; then
	sh ${TOOL_CONVERT_METH_CALLS_MOABS} ${NP_METHCALLS_CH} 'CH' > ${FILENAME_METH_CALL_CH_TMP} & procIDMethConvCH=$!;
	sh ${TOOL_CONVERT_METH_CALLS_MOABS} ${NP_METHCALLS_CG} 'CG' > ${FILENAME_METH_CALL_CG_TMP} & procIDMethConvCG=$!;
elif [[ ${METH_CALLS_CONVERTER} == 'none' ]]; then
	cat ${NP_METHCALLS_CH} | grep "CH" > ${FILENAME_METH_CALL_CH_TMP} & procIDMethConvCH=$!;
	cat ${NP_METHCALLS_CG} | grep "CG" > ${FILENAME_METH_CALL_CG_TMP} & procIDMethConvCG=$!;
else
    throw 50 "Unknown methylation call converter: METH_CALLS_CONVERTER='$METH_CALLS_CONVERTER'"
fi


# Apply "normal" B-Caller
if [ ${IS_TAGMENTATION-false} ]; then
	${PYTHON_BINARY} ${TOOL_METHYLATION_CALLING_SCRIPT} ${METH_CALL_PARAMETERS} \
		-r ${chr} \
		-m ${FILENAME_METH_CALL_METRICS} \
		${CYTOSINE_POSITIONS_INDEX} \
		${FILENAME_MERGED_BAM} - | \
	tee ${NP_METHCALLS_CG}  > ${NP_METHCALLS_CH} & procID_BCALLING=$!
# Apply Tagmentation specific B-Caller
else
	${PYTHON_BINARY} ${TOOL_METHYLATION_CALLING_SCRIPT_TAGMENTATION} ${METH_CALL_PARAMETERS} \
		-r ${chr} \
		-m ${FILENAME_METH_CALL_METRICS} \
		${CYTOSINE_POSITIONS_INDEX} \
		${FILENAME_MERGED_BAM} - | \
	tee ${NP_METHCALLS_CG}  > ${NP_METHCALLS_CH} & procID_BCALLING=$!
fi;

# Error Checking
wait $procIDMethConvCG || throw 14 "Error from MOABS meth call conversion (CG)"
wait $procIDMethConvCH || throw 15 "Error from MOABS meth call conversion (CH)"
wait $procID_BCALLING || throw 16 "Error in BCalling"

########################################
# Compress and Tabix methylation calls #
########################################
if [[ ${METH_CALLS_CONVERTER} == 'moabs' ]]; then
	( mv ${FILENAME_METH_CALL_CG_TMP} ${FILENAME_METH_CALL_CG} && \
		bgzip ${FILENAME_METH_CALL_CG} && \
		tabix -p bed ${FILENAME_METH_CALL_CG}.gz ) & procIDTabixCG=$!
	( mv ${FILENAME_METH_CALL_CH_TMP} ${FILENAME_METH_CALL_CH} && \
		bgzip ${FILENAME_METH_CALL_CH} && \
		tabix -p bed ${FILENAME_METH_CALL_CH}.gz ) & procIDTabixCH=$!
elif [[ ${METH_CALLS_CONVERTER} == 'none' ]]; then
	( mv ${FILENAME_METH_CALL_CG_TMP} ${FILENAME_METH_CALL_CG} && \
		bgzip ${FILENAME_METH_CALL_CG} && \
		tabix -s 1 -b 2 ${FILENAME_METH_CALL_CG}.gz ) & procIDTabixCG=$!
	( mv ${FILENAME_METH_CALL_CH_TMP} ${FILENAME_METH_CALL_CH} && \
		bgzip ${FILENAME_METH_CALL_CH} && \
		tabix -s 1 -b 2 ${FILENAME_METH_CALL_CH}.gz ) & procIDTabixCH=$!
fi;

# Error checking
wait ${procIDTabixCG} || throw 17 "Error from tabix indexing CG"
wait ${procIDTabixCH} || throw 18 "Error from tabix indexing CH"

touch ${FILENAME_METH_CALLS_CHECKPOINT}

exit 0
