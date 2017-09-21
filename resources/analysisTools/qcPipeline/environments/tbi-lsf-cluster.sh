#!/usr/bin/env bash

source "$TOOL_WORKFLOW_LIB"

module load R/3.0.0
module load bwa/"${BWA_VERSION:-0.7.8}"      # select version!
module load fastqc/1.11.3
module load perl/5.20.2
module load picard/1.125
module load python/2.7.9
module load samtools/0.1.19
module load htslib/0.2.5
module load bedtools/2.16.2

if markWithBiobamba; then
    module load libmaus/0.0.131   # used to be 0.0.130
    module load biobambam/0.0.148
fi

if markWithSambamba; then
    module load sambamba/"${SAMBAMBA_VERSION:-0.5.9}"
fi

if markWithPicard; then
    module load picard/1.125
fi

if [[ "WGBS" ]]; then
    module load moabs/1.3.0
fi

