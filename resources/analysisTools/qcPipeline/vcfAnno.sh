#!/bin/bash

set -o pipefail
source "$CONFIG_FILE"
source "$TOOL_BASH_LIB"
set -x

# Estimate gender of patient from X and Y coverage
if isControlSample "$SAMPLE"
then
	tmp_sex_file="${FILENAME_SEX}_tmp"

	${RSCRIPT_BINARY} "$TOOL_ESTIMATE_SEX" \
		 --file_size "$CHROMOSOME_LENGTH_FILE" \
		 --cnv_file "$FILENAME_COV_WINDOWS_1KB" \
		 --min_Y_ratio "$min_Y_ratio" \
		 --min_X_ratio "$min_X_ratio" \
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

# Reformat original coverage file so it can be annotated with annotate_vcf.pl.
A_FILE="${FILENAME_COV_WINDOWS_1KB}.end.txt.tmp"
cat "$FILENAME_COV_WINDOWS_1KB" | awk '{print $1,$2,$3+999,$3}' | sed 's/ /\t/g' | sed 's/^/chr/' |  sed '1i\#chr\tpos\tend\tcoverage' > "$A_FILE"

O_FILE="$FILENAME_COV_WINDOWS_1KB_ANNO"
tmp_out="${O_FILE}_tmp"

${PERL_BINARY} "$TOOL_ANNOTATE_CNV_VCF" \
      -a "$A_FILE" \
      --aFileType=custom \
      --aChromColumn chr \
      --aPosColumn pos \
      --aEndColumn end \
      -b "$MAPPABILITY_FILE" \
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


 

