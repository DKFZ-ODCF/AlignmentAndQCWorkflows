#!/usr/bin/sh

set -o pipefail
source "$TOOL_WORKFLOW_LIB"
set -x

tmp_corrected_windowfile="${FILENAME_GC_CORRECTED_WINDOWS}.tmp"
tmp_corrected_table_slim="${FILENAME_GC_CORRECTED_QUALITY}.tmp"

${RSCRIPT_BINARY} "$TOOL_CORRECT_GC_BIAS_R" \
	--windowFile	"$FILENAME_COV_WINDOWS_WG" \
	--timefile	"$REPLICATION_TIME_FILE_ALN" \
	--chrLengthFile	"$CHROMOSOME_LENGTH_FILE_ALN" \
	--pid		"$PID" \
	--sample	$(sampleType "$SAMPLE") \
	--outfile	"$tmp_corrected_windowfile" \
	--corPlot	"$FILENAME_GC_CORRECT_PLOT" \
	--qcTab		"$tmp_corrected_table_slim" \
	--gcFile	"$GC_CONTENT_FILE_ALN" \
	--outDir	"$aceseqOutputDirectory" \
	--lowess_f	"$LOWESS_F_ALN" \
	--scaleFactor	"$SCALE_FACTOR_ALN" \
	--coverageYlims "$COVERAGEPLOT_YLIMS_ALN"


if [[ $? != 0 ]]
then
	echo "Something went wrong during GC correction. Program had non-zero exit status, exiting pipeline...\n\n"
	exit 2
fi	

mv "$tmp_corrected_windowfile" "$FILENAME_GC_CORRECTED_WINDOWS"

tmp_qc_value_file="$FILENAME_QC_GC_CORRECTION_JSON.tmp"

${PYTHON_BINARY} "$TOOL_CONVERT_TO_JSON" \
	--file 	"$tmp_corrected_table_slim" \
	--id 	"$pid" \
	--out 	"$tmp_qc_value_file" \
	--key	"$GC_bias_json_key_ALN"

if [[ $? != 0 ]]
then
	echo "Something went wrong during json file conversion. Program had non-zero exit status, exiting pipeline...\n\n"
	exit 2
fi	

mv "$tmp_qc_value_file" "$FILENAME_QC_GC_CORRECTION_JSON"
mv "$tmp_corrected_table_slim" "$FILENAME_GC_CORRECTED_QUALITY"
