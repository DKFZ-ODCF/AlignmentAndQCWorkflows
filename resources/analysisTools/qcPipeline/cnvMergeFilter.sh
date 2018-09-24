#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

set -o pipefail
set -x

${PYTHON_BINARY} "${TOOL_MERGE_FILTER_CNV}" \
            --inputfile    "$FILENAME_COV_WINDOWS_ANNO" \
            --output       "$FILENAME_COV_WINDOWS_WG" \
            --coverage     "$cnv_min_coverage_ALN" \
            --mappability  "$mapping_quality_ALN" \
            --NoOfWindows  "$min_windows_ALN"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code in the cnv merge and filter process;" 
	exit 2
fi

