#!/bin/bash

set -o pipefail
source "$CONFIG_FILE"
set -x

${PYTHON_BINARY} "${TOOL_MERGE_FILTER_CNV}" \
            --inputfile    "$FILENAME_COV_WINDOWS_1KB_ANNO" \
            --output       "$FILENAME_COV_WINDOWS_WG" \
	        --coverage     "$cnv_min_coverage" \
            --mappability  "$mapping_quality" \
            --NoOfWindows  "$min_windows"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the cnv merge and filter process;" 
	exit 2
fi

