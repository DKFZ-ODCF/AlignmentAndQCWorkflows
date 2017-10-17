#!/usr/bin/env bash

TESTFILE=${DIR_SNPCOMP}/snpout_${LANE}_${PBS_ARRAYID}.txt

OUTFILE=${DIR_SNPCOMP}/${LANE}_vs_controlSNP_${snpCnt0}x${snpCnt1}_#CHROM#.txt
WARNING_FILE=${DIR_SNPCOMP}/${LANE}_vs_controlSNP_${snpCnt0}x${snpCnt1}_WARNING

echo "$OUTFILE, $WARNING_FILE" >> $TESTFILE

# convert to chromosome
#CHROM="chr$CHROM_INDEX"
CHROM=chr${PBS_ARRAYID}
[ $CHROM = chr23 ] && CHROM=chrX
[ $CHROM = chr24 ] && CHROM=chrY

#Affy SNP array coordinates (443011 lines) for target pileup    # made from /icgc/lsdf/mb/analysis/wangq/AffySNPArray/GPL6804-20642.hg19.sorted.bed with /home/hutter/Perlprog/tools/SNP_bed2vcf.pl
SNP_REFERENCE_ANNOTATIONS=${SNP_REFERENCE_ANNOTATIONS/\#CHROM#/$CHROM} #/icgc/lsdf/mb/analysis/Annotation/hg19/Affy5/${CHROM}_AFFY.vcf
OUTFILE=${OUTFILE/\#CHROM#/$CHROM}

${SAMTOOLS_BINARY} mpileup -RI -q 1 -l ${SNP_REFERENCE_ANNOTATIONS} -r $CHROM -f ${SNP_REFERENCE} ${f0} ${f1} | ${PYTHON_BINARY} ${TOOL_INSILICO_GENOTYPER_SCRIPT} - $OUTFILE ${SNP_MINCOVERAGE} ${SNP_MAXCOVERAGE}

#TODO Email mit Warning Datei
echo "done"  >> $TESTFILE
echo "" >> $TESTFILE
