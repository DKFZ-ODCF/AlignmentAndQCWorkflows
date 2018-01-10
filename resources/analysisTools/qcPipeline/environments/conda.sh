#!/usr/bin/env bash

# Copyright (c) 2017
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
# Load a Conda environment.

createCleanCondaEnvironment () {
    unset LD_LIBRARY_PATH
    unset MODULE_PATH
    unset PKG_CONFIG_PATH
    unset PYTHONHOME
    unset PYTHONPATH
    unset PYTHONSTARTUP
    unset PYTHON_LIB
    unset PBS_SHARED_BIN
    unset PERL5LIB
    unset PERL_LOCAL_LIB_ROOT
    unset PERL_MB_OPT
    unset PERL_MM_OPT
    export PATH="$HOME/miniconda3/bin:/opt/torque/bin:/usr/lib64/mpi/gcc/openmpi/bin:/opt/maui/bin:/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/X11R6/bin:/usr/lib/mit/bin"
}

createCleanCondaEnvironment


source activate "${condaEnvironmentName:?No Conda environment name defined. Please set 'condaEnvironmentName'.}" \
    || (echo "Could not load Conda environment '$condaEnvironmentName'" && exit 100)

export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export RSCRIPT_BINARY=Rscript
export JAVA_BINARY=java
export FASTQC_BINARY=fastqc
export PERL_BINARY=perl
export PYTHON_BINARY=python
export SAMTOOLS_BINARY=samtools
export BCFTOOLS_BINARY=bcftools
export INTERSECTBED_BINARY=intersectBed
export COVERAGEBED_BINARY=coverageBed
export FASTAFROMBED_BINARY=fastaFromBed
export BAMSORT_BINARY=bamsort
export BAMMARKDUPLICATES_BINARY=bammarkduplicates
export PICARD_BINARY=picard.sh
export VCFTOOLS_SORT_BINARY=vcf-sort
export SAMBAMBA_BINARY=sambamba
export SAMBAMBA_FLAGSTATS_BINARY=sambamba
export SAMBAMBA_MARKDUP_BINARY=sambamba
export MBUFFER_BINARY=mbuffer
export CHECKSUM_BINARY=md5sum
export BWA_BINARY=bwa
