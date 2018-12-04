#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#
# Load a Conda environment.

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
export TRIMMOMATIC_BINARY=trimmomatic
