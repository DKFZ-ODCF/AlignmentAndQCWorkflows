#!/usr/bin/env bash

set -xv


## From http://unix.stackexchange.com/questions/26676/how-to-check-if-a-shell-is-login-interactive-batch
shellIsInteractive () {
    case $- in
        *i*) echo "true";;
        *)   echo "false";;
    esac
}
export -f shellIsInteractive


## funname () ( set +exv; ...; ) may be better to get rid of too much output (mind the (...) subshell) but the exit won't work anymore.
## Maybe set -E + trap "bla" ERR would work? http://fvue.nl/wiki/Bash:_Error_handling#Exit_on_error
printStackTrace () {
    frameNumber=0
    while caller $frameNumber ;do
      ((frameNumber++))
    done
}
export -f printStackTrace


errout () {
    local exitCode="$1"
    local message="$2"
    env printf "Error(%d): %s\n" "$exitCode" "$message" >> /dev/stderr
}
export -f errout


## This is to effectively debug on the command line. The exit is only called, in non-interactive sessions.
## You can either put 'exitIfNonInteractive $code; return $?' at the end of functions, or you put
## 'exitHere $code || return $?' in the middle of functions to end the control flow in the function and
## return to the calling function.
exitIfNonInteractive () {
    local exitValue="$1"
    if [[ $(shellIsInteractive) == false ]]; then
      exit "$exitValue"
    else
      echo "In a non-interactive session, I would now do 'exit $exitValue'" >> /dev/stderr
      return "$exitValue"
    fi
}
export -f exitIfNonInteractive

## throw [code [msg]]
## Write message (Unspecified error) to STDERR and exit with code (default 1)
throw () {
  local lastCommandsExitCode=$?
  local exitCode="${1-$UNSPECIFIED_ERROR_CODE}"
  local msg="${2-$UNSPECIFIED_ERROR_MSG}"
  if [[ $lastCommandsExitCode -ne 0 ]]; then
    msg="$msg (last exit code: $lastCommandsExitCode)"
  fi
  errout "$exitCode" "$msg"
  printStackTrace
  exitIfNonInteractive "$exitCode" || return $?
}
export -f throw

# Load/Unload a module with $name and using the version given by the $versionVariable.
# If the version is not given, take the name, put it into upper case and append _VERSION.
versionVariable () {
    local name="${1:-No binary name given}"
    local versionVariable="$2"
    if [[ -z "$versionVariable" ]]; then
        local ucVersionVariable
        ucVersionVariable=$(echo "$name" | tr '[a-z]' '[A-Z]')
        echo "${ucVersionVariable}_VERSION"
    else
        echo "$versionVariable"
    fi
}
export -f versionVariable

moduleLoad() {
    local name="${1:?No module name given}"
    local versionVariable
    versionVariable=$(versionVariable "$name" $2)
    if [[ -z "${!versionVariable}" ]]; then
        throw 200 "$versionVariable is not set" > /dev/stderr
    fi
    module load "${name}/${!versionVariable}" || exit 200
}
export -f moduleLoad

moduleUnload() {
    local name="${1:?No module name given}"
    local versionVariable=$(versionVariable "$name" $2)
    if [[ -z "${!versionVariable}" ]]; then
        throw 200 "versionVariable for $name is not set" > /dev/stderr
    fi
    module unload "${name}/${!versionVariable}" || exit 200
}
export -f moduleUnload

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
export -f sambamba_sort_view
export SAMBAMBA_BINARY=sambamba

# The sambamba version used only for flagstats. For the flagstats sambamba 0.4.6 used is equivalent to samtools 0.1.19 flagstats. Newer versions
# use the new way of counting in samtools (accounting for supplementary reads).
# Warning: Currently bwaMemSortSlim uses sambamba flagstats, while mergeAndMarkOrRemoveSlim uses samtools flagstats.
sambamba_flagstat() {
    moduleLoad sambamba SAMBAMBA_FLAGSTATS_VERSION
    sambamba "$@"
    moduleUnload sambamba SAMBAMBA_FLAGSTATS_VERSION
}
export -f sambamba_flagstat
export SAMBAMBA_FLAGSTATS_BINARY=sambamba_flagstat

# The sambamba version used only for duplication marking and merging. Use the bash function here!
# Should be changeable independently also for performance and stability reasons.
sambamba_markdup() {
    moduleLoad sambamba SAMBAMBA_MARKDUP_VERSION
    sambamba "$@"
    moduleUnload sambamba SAMBAMBA_MARKDUP_VERSION
}
export -f sambamba_markdup
export SAMBAMBA_MARKDUP_BINARY=sambamba_markdup


if [[ "$WORKFLOW_ID" == "bisulfiteWorkflow" ]]; then
    ## For bisulfite alignment, we suffix the the value of BINARY_VERSION by '-bisulfite', because that's the name in LSF cluster.
    export BWA_VERSION="${BWA_VERSION:?BWA_VERSION is not set}-bisulfite"
    moduleLoad bwa
    export BWA_BINARY=bwa

    moduleLoad
elif [[ "$WORKFLOW_ID" == "qcPipeline" || "$WORKFLOW_ID" == "exomePipeline" ]]; then
    if [[ "${useAcceleratedHardware:-false}" == false ]]; then
        moduleLoad bwa
        export BWA_BINARY=bwa
    elif [[ "${useAcceleratedHardware:-true}" == true ]]; then
        moduleLoad bwa-bb BWA_VERSION
        export BWA_ACCELERATED_BINARY=bwa-bb
    else
        throw 200 "Uninterpretable value for boolean 'useAcceleratedHardware': '$useAcceleratedHardware'"
    fi
else
    throw 200 "Unknown workflow ID '$WORKFLOW_ID'"
fi



#    if [[ ACCELERATED OR FPGA && FPGA-Node? ]]; then
#        moduleLoad bwa-bb BWA_VERSION
#        export BWA_ACCELERATED_BINARY=bwa-bb
#        export BWA_BINARY=bwa-bb
#    else
#        moduleLoad bwa
#        export BWA_BINARY=bwa
#    fi

# Unversioned binaries.
export MBUFFER_BINARY=mbuffer
export CHECKSUM_BINARY=md5sum
