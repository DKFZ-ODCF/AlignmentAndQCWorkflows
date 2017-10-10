#!/usr/bin/env bash

set -vex

# Load a module with $name and using the version given by the $versionVariable.
# If the version is not given, take the name, put it into upper case and append _VERSION
versionVariable () {
    local versionVariable="$1"
    if [[ -z "$versionVariable" ]]; then
        eche $(echo "$name" | tr [a-z] [A-Z])"_VERSION"
    else
        echo "$versionVariable"
    fi
}
moduleLoad() {
    local name="${1:?No module name given}"
    local versionVariable=$(versionVariable "$2")
    module load "${name}/${!versionVariable}" || exit 200
}
moduleUnload() {
    local name="${1:?No module name given}"
    local versionVariable=$(versionVariable "$2")
    module unload "${name}/${!versionVariable}" || exit 200
}

moduleLoad htslib
export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix

moduleLoad R
export RSCRIPT_BINARY=Rscript

moduleLoad java
export JAVA_BINARY=java

moduleLoad fastqc
export FASTQC_BINARY=fastqc

moduleLoad perl
export PERL_BINARY=perl

moduleLoad python
export PYTHON_BINARY=python

moduleLoad pypy
export PYPY_BINARY=pypy-c

moduleLoad samtools
export SAMTOOLS_BINARY=samtools
export BCFTOOLS_BINARY=bcftools

moduleLoad bedtools
export INTERSECTBED_BINARY=intersectBed
export COVERAGEBED_BINARY=coverageBed
export FASTAFROMBED_BINARY=fastaFromBed

moduleLoad libmaus
moduleLoad biobambam
export BAMSORT_BINARY=bamsort
export BAMMARKDUPLICATES_BINARY=bammarkduplicates

moduleLoad picard
export PICARD_BINARY=picard.sh

moduleLoad vcftools
export VCFTOOLS_SORT_BINARY=vcf-sort

# There are different sambamba versions used for different tasks. The reason has something to do with performance and differences in the versions.
# We define functions here that take the same parameters as the original sambamba, but apply them to the appropriate version by first loading
# the correct version and after the call unloading the version. The _BINARY variables are set to the functions.

# The sambamba version used for sorting, viewing. Note that v0.5.9 is segfaulting on Convey during view or sort.
sambamba_sort_view() {
    moduleLoad sambamba
    sambamba "$@"
    moduleUnload sambamba
}
export SAMBAMBA_BINARY=sambamba

# The sambamba version used only for flagstats. For the flagstats sambamba 0.4.6 used is equivalent to samtools 0.1.19 flagstats. Newer versions
# use the new way of counting in samtools (accounting for supplementary reads).
# Warning: Currently bwaMemSortSlim uses sambamba flagstats, while mergeAndMarkOrRemoveSlim uses samtools flagstats.
sambamba_flagstat() {
    moduleLoad sambamba SAMBAMBA_FLAGSTAT_VERSION
    sambamba "$@"
    moduleUnload sambamba SAMBAMBA_FLAGSTAT_VERSION
}
export SAMBAMBA_FLAGSTATS_BINARY=sambamba_flagstat

# The sambamba version used only for duplication marking and merging. Use the bash function here!
# Should be changeable independently also for performance and stability reasons.
sambamba_markdup() {
    moduleLoad sambamba SAMBAMBA_MARKDUP_VERSION
    sambamba "$@"
    moduleUnload sambamba SAMBAMBA_MARKDUP_VERSION
}
export SAMBAMBA_MARKDUP_BINARY=sambamba_markdup

# Dependent on whether alignments are done on the FPGA or whether WGBS data are processed different BWA versions need to be loaded.
if [[ "$WORKFLOW_ID" == "bisulfiteCoreAnalysis" ]]; then
    # This is actually only a lightly patched version that does not check for read numbers in identifiers.
    moduleLoad bwa BWA_BISULPHITE_VERSION
    export BWA_BINARY=bwa
else
    if [[ ACCELERATED OR FPGA && FPGA-Node? ]]; then
        moduleLoad bwa-bb BWA_VERSION
        export BWA_ACCELERATED_BINARY=bwa-bb
        export BWA_BINARY=bwa-bb
    else
        moduleLoad bwa
        export BWA_BINARY=bwa
    fi

fi

# Unversioned binaries.
export MBUFFER_BINARY=mbuffer
export CHECKSUM_BINARY=md5sum
