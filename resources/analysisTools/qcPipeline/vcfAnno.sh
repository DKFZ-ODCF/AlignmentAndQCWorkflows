#!/bin/bash

set -o pipefail
source "$TOOL_WORKFLOW_LIB"
set -x

# Reformat original coverage file so it can be annotated with annotate_vcf.pl.
A_FILE="${FILENAME_COV_WINDOWS}.end.txt.tmp"
cat "$FILENAME_COV_WINDOWS" | awk '{print $1,$2,$2+999,$3}' | sed 's/ /\t/g' | sed 's/^/chr/' |  sed '1i\#chr\tpos\tend\tcoverage' > "$A_FILE"

# Estimate gender of patient from X and Y coverage
if isControlSample "$SAMPLE"
then
	tmp_sex_file="${FILENAME_SEX}_tmp"

	${RSCRIPT_BINARY} "$TOOL_ESTIMATE_SEX" \
		 --file_size "$CHROMOSOME_LENGTH_FILE_ALN" \
		 --cnv_file "$A_FILE" \
		 --min_Y_ratio "$min_Y_ratio_ALN" \
		 --min_X_ratio "$min_X_ratio_ALN" \
		 --file_out "$tmp_sex_file"


	if [[ "$?" != 0 ]]
	then
		echo "There was a non-zero exit code in the gender estimation" 
		exit 2
	fi

	mv "$tmp_sex_file" "$FILENAME_SEX"
else
    echo "Gender-determination unsafe for tumor sample '$SAMPLE'." > "$FILENAME_SEX"
fi

O_FILE="$FILENAME_COV_WINDOWS_ANNO"
tmp_out="${O_FILE}_tmp"

${PERL_BINARY} "$TOOL_ANNOTATE_CNV_VCF" \
      -a "$A_FILE" \
      --aFileType=custom \
      --aChromColumn chr \
      --aPosColumn pos \
      --aEndColumn end \
      -b "$MAPPABILITY_FILE_ALN" \
      --bFileType=bed \
      --reportBFeatCoord \
      --columnName map |
${PYTHON_BINARY} "$TOOL_ADD_MAPPABILITY" \
      -o "$tmp_out"

if [[ "$?" != 0 ]]
then
	echo "There was a non-zero exit code while adding the mappability values."
	exit 2
fi

mv "$tmp_out" "$O_FILE"
rm "$A_FILE"


 

