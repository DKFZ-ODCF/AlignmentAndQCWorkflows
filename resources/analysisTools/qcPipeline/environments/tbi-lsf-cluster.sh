#!/usr/bin/env bash

SAMBAMBA_VERSION=0.5.9
SAMBAMBA_FLAGSTAT_VERSION=0.4.6
SAMBAMBA_MARKDUP_VERSION=0.5.9
module load R/3.0.0
module load bwa/"${BWA_VERSION:-0.7.8}"      # select version!
module load java/1.8.0_131
module load fastqc/1.11.3
module load perl/5.20.2
module load python/2.7.9
module load samtools/0.1.19
module load htslib/0.2.5
module load bedtools/2.16.2

if markWithBiobambam; then
    module load libmaus/0.0.130
    module load biobambam/0.0.148
elif markWithPicard; then
    module load picard/1.125
else
    echo "Oops" > /dev/stderr
fi

if [[ "WGBS" ]]; then
    module load moabs/1.3.0
fi

# The sambamba version used for sorting, viewing. Note that v0.5.9 is segfaulting on convey during view or sort.
sambamba_sort_view() {
    module load "sambamba/$SAMBAMBA_VERSION"
    sambamba "$@"
}
export SAMBAMBA_BINARY=sambamba

# The sambamba version used only for making flagstats.
# Warning: Currently bwaMemSortSlim uses sambamba flagstats, while mergeAndMarkOrRemoveSlim uses samtools flagstats.
sambamba_flagstat() {
    module load "sambamba/$SAMBAMBA_FLAGSTAT_VERSION"
    sambamba "$@"
}
export SAMBAMBA_FLAGSTATS_BINARY=sambamba_flagstat

# The sambamba version used only for duplication marking and merging. Use the bash function here!
sambamba_markdup() {
    module load "sambamba/$SAMBAMBA_MARKDUP_VERSION"
    sambamba "$@"
}
export SAMBAMBA_MARKDUP_BINARY=sambamba_markdup

export PICARD_BINARY=picard.sh
export FASTQC_BINARY=fastqc
export SAMTOOLS_BINARY=samtools

# biobambam
export BAMSORT_BINARY=bamsort
export MARKDUPLICATES_BINARY=bammarkduplicates

export JAVA_BINARY=java

# htslib/0.2.5
export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix

export PERL_BINARY=perl
export PYTHON_BINARY=python
export PYPY_BINARY=pypy-c
export RSCRIPT_BINARY=Rscript
export MBUFFER_BINARY=mbuffer

#export FASTAFROMBED_BINARY=fastaFromBed-2.16.2
#export BCFTOOLS_BINARY=bcftools-0.1.19
#export VCFTOOLS_SORT_BINARY=vcf-sort-0.1.10

# WES only
export INTERSECTBED_BINARY=intersectBed
export COVERAGEBED_BINARY=coverageBed

if [[ WGBS ]]; then
    export BWA_BINARY=bwa-0.7.8-bisulfite
else
    export BWA_BINARY=bwa
    export BWA_ACCELERATED_BINARY=/opt/bb/bwa-0.7.8-r2.05/bin/bwa-bb
fi

