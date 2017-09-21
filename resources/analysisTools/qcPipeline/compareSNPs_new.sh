#!/bin/bash

#PBS -l walltime=2:00:00,nodes=1:ppn=1


TESTFILE=${DIR_SNPCOMP}/snpout_${LANE}.txt

# wie werden die Dateien genannt, und welche BAMs gehen ein?
# bei einzelnen lanes is die coverage niedriger und concordance auch, hier meist 95% und warnings
# dagegen 30x BAMs > 98% concordance (ausser bei LOH)
# daher prüfen, was es ist, und nur aussteigen, wenn > x Zeilen in der _WARNING Datei
OUTFILE=${DIR_SNPCOMP}/${LANE}_vs_controlSNP_${snpCnt0}x${snpCnt1}.txt
WARNING_FILE=${DIR_SNPCOMP}/${LANE}_vs_controlSNP_${snpCnt0}x${snpCnt1}_WARNING

echo "$OUTFILE, $WARNING_FILE" >> $TESTFILE

for chridx in ${CHROMOSOME_INDICES[@]}
do
	chrom=${CHR_PREFIX}$chridx
	AFFY_VCF=${SNP_REFERENCE_ANNOTATIONS/\#CHROM#/$chrom}
	if [ ! -e $AFFY_VCF ]
	then
		echo "$AFFY_VCF does not exist, exiting"
		exit 22
	fi
	${SAMTOOLS_BINARY} mpileup -RI -q 1 -l $AFFY_VCF -r $chrom -f ${REFERENCE_GENOME} ${f0} ${f1} | ${PYTHON_BINARY} ${TOOL_INSILICO_GENOTYPER_SCRIPT} - $OUTFILE ${SNP_MINCOVERAGE} ${SNP_MAXCOVERAGE}


done

#TODO Email mit Warning Datei
echo "done"  >> $TESTFILE
echo "" >> $TESTFILE
