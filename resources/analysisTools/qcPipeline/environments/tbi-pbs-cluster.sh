#!/usr/bin/env bash

# The sambamba version used for sorting, viewing. Note that v0.5.9 is segfaulting on convey during view or sort.
export SAMBAMBA_BINARY=sambamba_v0.4.6

# The sambamba version used only for making flagstats.
export SAMBAMBA_FLAGSTATS_BINARY=sambamba_v0.4.6
# Warning: Currently bwaMemSortSlim uses sambamba flagstats, while mergeAndMarkOrRemoveSlim uses samtools flagstats.

# The sambamba version used only for duplication marking and merging.
export SAMBAMBA_MARKDUP_BINARY=sambamba_v0.5.9

export PICARD_BINARY=picard-1.125.sh
export FASTQC_BINARY=fastqc-0.10.1
export SAMTOOLS_BINARY=samtools-0.1.19

# biobambam
export BAMSORT_BINARY=bamsort-0.0.148
export BAMMARKDUPLICATES_BINARY=bammarkduplicates-0.0.148

export JAVA_BINARY=java8

# htslib/0.2.5
export BGZIP_BINARY=bgzip-tabix-0.2.5
export TABIX_BINARY=tabix-0.2.5

export PERL_BINARY=perl-5.20.2
export PYTHON_BINARY=python-2.7.9
export PYPY_BINARY=pypy-c
export RSCRIPT_BINARY=Rscript-3.0.0
export MBUFFER_BINARY=mbuffer

#export FASTAFROMBED_BINARY=fastaFromBed-2.16.2
#export BCFTOOLS_BINARY=bcftools-0.1.19
#export VCFTOOLS_SORT_BINARY=vcf-sort-0.1.10

# WES only
export INTERSECTBED_BINARY=intersectBed-2.16.2
export COVERAGEBED_BINARY=coverageBed-2.16.2

if [[ WGBS ]]; then
    export BWA_BINARY=bwa-0.7.8-bisulfite
else
    export BWA_BINARY=bwa-0.7.8
    export BWA_ACCELERATED_BINARY=/opt/bb/bwa-0.7.8-r2.05/bin/bwa-bb
fi



