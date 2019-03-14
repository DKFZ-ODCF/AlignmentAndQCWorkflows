#
# Copyright (c) 2019 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#
##############################################################################
## Domain-specific code
##############################################################################
##
## Include into your code with: source "$TOOL_WORKFLOW_LIB"
##

source "$TOOL_BASH_LIB"

WORKFLOWLIB___ERREXIT=$(if [[ $SHELLOPTS =~ "errexit" ]]; then echo "errexit"; fi)
WORKFLOWLIB___SHELL_OPTIONS=$(set +o)
set +o verbose
set +o xtrace

markWithPicard () {
    [[ "$markDuplicatesVariant" == "picard" || ("$markDuplicatesVariant" == "" && "$useBioBamBamMarkDuplicates" == "false") ]]
    return $?
}


markWithSambamba () {
    [[ "$markDuplicatesVariant" == "sambamba" ]]
    return $?
}


markWithBiobambam () {
    [[ "$markDuplicatesVariant" == "biobambam" || ("$markDuplicatesVariant" == "" && "$useBioBamBamMarkDuplicates" == "true") ]]
    return $?
}


mbuf () {
    local bufferSize="$1"
    shift
    assertNonEmpty "$bufferSize" "No buffer size defined for mbuf()" || return $?
    "$MBUFFER_BINARY" -m "$bufferSize" -q -l /dev/null "$@"
}


## The return the directory, to which big temporary files should be written, e.g. for sorting.
getBigScratchDirectory () {
    local suggestedLocation="${1:-}"
    local scratchDir
    if [[ "${useRoddyScratchAsBigFileScratch:-false}" == true ]]; then
        scratchDir="${RODDY_SCRATCH}"
    elif [[ -z $suggestedLocation ]]; then
        scratchDir="$outputAnalysisBaseDirectory/tmp"
    else
        scratchDir="$suggestedLocation"
    fi
    mkdir -p "$scratchDir" || throw "$NOT_WRITABLE_CODE" "$NOT_WRITABLE_MSG: '$scratchDir'"
    echo "$scratchDir"
}



analysisType () {
    if [[  "${runExomeAnalysis-false}" = "true" ]]; then
        echo "exome"
    else
        echo "genome"
    fi
}


md5File () {
   local inputFile="${1-/dev/stdin}"
   local outputFile="${2-/dev/stdout}"
   assertNonEmpty "$inputFile"  "inputFile not defined" || return $?
   assertNonEmpty "$outputFile" "outputFile not defined" || return $?
   cat "$inputFile" \
        | ${CHECKSUM_BINARY} \
        | cut -d ' ' -f 1 \
        > "$outputFile"
}


samtoolsIndex () {
   local inputFile="${1-/dev/stdin}"
   local outputFile="${2-/dev/stdout}"
   assertNonEmpty "$inputFile"  "inputFile not defined" || return $?
   assertNonEmpty "$outputFile" "outputFile not defined" || return $?
   "$SAMTOOLS_BINARY" index "$inputFile" "$outputFile"
}


fakeDupMarkMetrics () {
   local inputFile="${1-/dev/stdin}"
   local outputFile="${2-/dev/stdout}"
   assertNonEmpty "$inputFile"  "inputFile not defined" || return $?
   assertNonEmpty "$outputFile" "outputFile not defined" || return $?
   "$PERL_BINARY" "$TOOL_FAKE_DUPMARK_METRICS" "$inputFile" "${SAMPLE}_${pid}" \
        > "$outputFile"
}


toIEqualsList () {
    declare -a inputFiles=($@)
    for inFile in ${inputFiles[@]}; do
	    echo -n "I=$inFile "
    done
    echo
}


# Stolen from http://stackoverflow.com/questions/3685970/check-if-an-array-contains-a-value
arrayContains () {
    local ELEMENT="${1}"
    assertNonEmpty "$ELEMENT" "arrayContains called without parameters" || return $?
    local DELIM=","
    printf "${DELIM}%s${DELIM}" "${@:2}" | grep -q "${DELIM}${ELEMENT}${DELIM}"
}

matchPrefixInArray () {
    local ELEMENT="${1}"
    assertNonEmpty "$ELEMENT" "matchPrefixInArray called without parameters" || return $?
    local DELIM=" "
    printf "${DELIM}%s${DELIM}" "${@:2}" | grep -q -P "${DELIM}${ELEMENT}[^${DELIM}]*${DELIM}"
}

isControlSample () {
    assertNonEmpty "$1" "isControlSample expects sample type name as single parameter" || return $?
    declare -a prefixes="${possibleControlSampleNamePrefixes[@]}"
    matchPrefixInArray "$1" "${prefixes[@]}"
}

isTumorSample () {
    assertNonEmpty "$1" "isTumorSample expects sample type name as single parameter" || return $?
    declare -a prefixes="${possibleTumorSampleNamePrefixes[@]}"
    matchPrefixInArray "$1" "${prefixes[@]}"
}

sampleType () {
    assertNonEmpty "$1" "sampleType expects sample type name as single parameter" || return $?
    if isControlSample "$1" && isTumorSample "$1"; then
        throw_illegal_argument "Sample '$1' cannot be control and tumor at the same time"
    elif isControlSample "$1"; then
        echo "control"
    elif isTumorSample "$1"; then
        echo "tumor"
    else
        throw_illegal_argument "'$1' is neither control nor tumor"
    fi
}

runFingerprinting () {
    local bamFile="${1:?No input BAM file given}"
    local fingerPrintsFile="${2:?No output fingerprints file given}"
    if [[ "${runFingerprinting:-false}" == true && ${useAcceleratedHardware:-false} != true ]]; then
        "${PYTHON_BINARY}" "${TOOL_FINGERPRINT}" "${fingerprintingSitesFile}" "${bamFile}" > "${fingerPrintsFile}.tmp" || throw 43 "Fingerprinting failed"
        mv "${fingerPrintsFile}.tmp" "${fingerPrintsFile}" || throw 39 "Could not move file"
    fi
}


# Remove a single directory, owned by the current user, recursively. Certain file names are forbidden and it is checked
# that only a single file or directory (including contained files) will be removed.
saferRemoveSingleDirRecursively () {
    local file="${1:?No file given}"
    if [[ "${file}" == "" || "${file}" == "." || "${file}" == "/" || "${file}" == "*" ]]; then
        throw 1 "Trying to recursively remove with forbidden file name: '${file}'"
    fi
    declare -a owner=( $(stat -c "%U" "${file}") )
    if [[ "${#owner[@]}" -gt 1 ]]; then
        throw 1 "Trying to remove multiple files with file pattern: '${file}'"
    fi
    if [[ "${owner}" != $(whoami) ]]; then
        throw 1 "${file} is owned by ${owner}, so it won't be deleted by $(whoami)"
    fi
    rm -rf "${file}"
}

removeRoddyBigScratch () {
    if [[ "${RODDY_BIG_SCRATCH}" && "${RODDY_BIG_SCRATCH}" != "${RODDY_SCRATCH}" ]]; then  # $RODDY_SCRATCH is also deleted by the wrapper.
        saferRemoveSingleDirRecursively "${RODDY_BIG_SCRATCH}" # Clean-up big-file scratch directory. Only called if no error in wait or mv before.
    fi
}

checkBamIsComplete () {
    local bamFile="${1:?No BAM file given}"
    local result      # AFAIR the assigment and declaration together did not work; presumably because of a Bash bug!
    result=$("$TOOL_BAM_IS_COMPLETE" "$bamFile")
    if [[ $? ]]; then
        echo "BAM is terminated! $bamFile" >> /dev/stderr
    else
        throw 40 "BAM is not terminated! $bamFile"
    fi
}

throwIfFileInaccessible() {
    local filename="${1:?No file given}"
    if [[ ! -r "$filename" ]]; then
        throw 100 "Cannot access file: $filename"
    fi
}

throwIfMatching() {
    local filename="${1:?No BWA STDERR output file given}"
    local pattern="${2:?No pattern to match given}"
    local errorCode="${3:?No error code to report given}"
    if [[ -r "$filename" ]]; then
        throw 37 "Cannot access '$filename'"            # Because of the `|| true` below let's check the basics first.
    fi
    errorMsg=$(grep -P "$pattern" "$filename" || true)  # The `|| true` copes with grep's silly exit 1 upon no-match.
    if [[ ! -z "$errorMsg" ]]; then
        throw 36 "Found error in logfile '$filename': [$errorMsg]"
    fi
}

# Note: Checks are not collected. To get as precise information as possible the tests are sorted from more specific (in terms of diagnostic value)
#       to less specific.
checkBwaOutput() {
    local bamFile="${1:?No BAM file given}"
    local bwaOutput="${2:?No BWA STDERR output file given}"
    local sortLog="${3:?No sort log given}"
    local errorCodeFile="${4:?No error code file given}"

    throwIfFileInaccessible "$bamFile"
    throwIfFileInaccessible "$bwaOutput"
    throwIfFileInaccessible "$errorCodeFile"

    # Check for segfault messages
    throwIfMatching "$bwaOutput" " fault" 31

    # Barbara Aug 10 2015: I can't remember what bwa aln and sampe reported as "error".
    # bluebee bwa has "error_count" in bwa-0.7.8-r2.05; and new in bwa-0.7.8-r2.06: "WARNING:top_bs_ke_be_hw: dummy be execution, only setting error."
    # these are not errors that would lead to fail, in contrast to "ERROR: Bus error"
    throwIfMatching <(cat "$bwaOutput" | grep -v "error_count" | grep -v "dummy be execution") "error" 36

    # An assumption of the workflow is violated.
    throwIfMatching "$bwaOutput" "file has fewer sequences." 41

    # Check for BWA abortion.
    throwIfMatching "$bwaOutput" "Abort. Sorry." 37

    # samtools sort may complain about truncated temp files and for each line outputs
    # the error message. This happens when the same files are written to at the same time,
    # See http://sourceforge.net/p/samtools/mailman/samtools-help/thread/BAA90EF6FE3B4D45A7B2F6E0EC5A8366DA3AB5@USTLMLLYC102.rf.lilly.com/
    # This happens when the scheduler puts the same job on 2 nodes bc. The prefix for samtools-0.1.19 -o $prefix is constructed using the job ID.
    if [[ ! -z $sortLog ]] && [[ -f $sortLog ]]; then
        throwIfMatching "$sortLog" "is truncated. Continue anyway." 38
    else
        echo "there is no samtools sort log file" >> /dev/stderr
    fi

    if [[ $(cat "$errorCodeFile") != "0" ]]; then
        throw 32 "There was a non-zero exit code in bwa pipe (w/ pipefail): $(cat '$errorCodeFile'); exiting..."
    fi

    if [[ "2048" -gt `stat -c %s $bamFile` ]]; then
        throw 33 "Output file is too small!"
    fi

    checkBamIsComplete "$bamFile"

    echo all OK
}


readGroupsInBam() {
    local bamFile="${1:?No bam file given}"
    local readGroups=`${SAMTOOLS_BINARY} view -H ${bamFile} | grep "^@RG"`
    if [[ -z "$readGroups" ]]; then
        throw 23 "could not detect single lane BAM files in $bamFile, stopping here"
    fi
    echo "$readGroups"
}

getMissingReadGroups() {
    local pairedBamSuffix="${1:?No paired-bam suffix given}"
    local sample="${2:?No sample name given}"
    local bamFile="${3:?No existing bam file given}"
    shift 3
    declare -a inputFiles=($@)

    local existingBamReadGroups=`readGroupsInBam "$bamFile"`
    ## Note: This does not test or even complain, if the BAM (e.g. due to manual manipulation) contains lanes that are
    ##       NOT among the INPUT_FILES. TODO Add at least a warning upon unknown lanes in BAM.
    readGroupsToMerge=`perl ${TOOL_PRINT_MISSING_READ_GROUPS} $(stringJoin ":" ${inputFiles[@]}) "$existingBamReadGroups" "$pairedBamSuffix" "$sample"`
    if [[ "$?" != 0 ]]; then
        throw 24 "something went wrong with the detection of merged files in ${EXISTING_BAM}, stopping here"
    fi
    echo "$readGroupsToMerge"
}


matchesShortChromosomeName() {
    local val="${1:?No value to match short chromosome pattern against}"
    if [[ $(matchesLongChromosomeName "$val") == "true" ]]; then
        echo "false"
    else
        echo "true"
    fi
}

matchesLongChromosomeName() {
    local val="${1:?No value to match long chromosome pattern against}"
    if [[ "${CHR_PREFIX:-}" != "" || "${CHR_SUFFIX:-}" != "" ]]; then
        if [[ "${val##${CHR_PREFIX:-}*${CHR_SUFFIX:-}}" == "" ]]; then
            echo "true"
        else
            echo "false"
        fi
    else
        echo "false"
    fi
}

chromosomeIndices() {
    if [[ ! -v ___chromosomeIndices ]]; then
        declare -g -a ___chromosomeIndices=( $(cut -f 1 "$CHROM_SIZES_FILE") )
    fi
    echo "${___chromosomeIndices[@]}"
}

shortChromosomeGroupSpec() {
    declare -a chromosomeIndices=( $(chromosomeIndices) )
    echo -n "${CHR_GROUP_NOT_MATCHING:-nonmatching}="
    declare -a shorts=()
    for chr in "${chromosomeIndices[@]}"; do
        if [[ $(matchesShortChromosomeName "$chr") == "true" ]]; then
            set +u
            shorts=("${shorts[@]}" "$chr")
        fi
    done
    echo $(set +u; stringJoin "," "${shorts[@]}")
}

longChromosomeGroupSpec() {
    declare -a chromosomeIndices=( $(chromosomeIndices) )
    echo -n "${CHR_GROUP_MATCHING:-matching}="
    declare -a longs=()
    for chr in "${chromosomeIndices[@]}"; do
        if [[ $(matchesLongChromosomeName "$chr") == "true" ]]; then
            set +u
            longs=("${longs[@]}" "$chr")
        fi
    done
    echo $(set +u; stringJoin "," "${longs[@]}")
}

groupLongAndShortChromosomeNames() {
    local genomeCoverageFile="${1:-/dev/stdin}"
    declare -a chromosomeIndices=( $(chromosomeIndices) )
    $PERL_BINARY $TOOL_GROUPED_GENOME_COVERAGES \
        "$genomeCoverageFile" \
        $(shortChromosomeGroupSpec) \
        $(longChromosomeGroupSpec)
}

#############################################################################
# Linear pipe management
#############################################################################

# Make a name for a read and a given pipe-name.
mkPairedPipeName() {
    local readNo="${1:?No 1/2 as read identifier}"
    local pipeName="${2:?No pipe basename}"
    if [[ $readNo -ne 1 && $readNo -ne 2 ]]; then
        throw 151 "Read number must be 1 or 2. Got '$readNo'"
    fi
    echo "r${readNo}_${pipeName}"
}

# Make and register a pair of named pipes from the pipe-name.
mkPipePairSource() {
    local pipeName="${1:?Named-pipe basename}"
    mkPipeSource $(mkPairedPipeName 1 "$pipeName")
    mkPipeSource $(mkPairedPipeName 2 "$pipeName")
}

# Get the path to the pipe named pipename for read 1 or 2.
getPairedPipeEndPath() {
    local readNo="${1:?No read-number}"
    local pipeName="${2:?Named-pipe basename}"
    getPipeEndPath $(mkPairedPipeName "$readNo" "$pipeName")
}

# Extend a pair of pipes by running a command with the interface "command inpipe1 inpipe2 outpipe1 outpipe2 @rest".
# The pipeline end is automatically progressed. Call like this:
#
#   extendPipePair $pipeName $tag -- commandName @restArgs
#
# The "--" is optional.
extendPipePair() {
    local pipeName="${1:?Named-pipe base path}"
    local tag="${2:?No processing step tag}"
    shift 2
    if [[ "$1" == "--" ]]; then
        shift
    fi
    local command="${1:?No command/function}"
    shift
    declare -a rest=("$@")

    local pipe1Name=$(mkPairedPipeName 1 "$pipeName")
    local pipe2Name=$(mkPairedPipeName 2 "$pipeName")

    local r1_inpipe=$(getPipeEndPath "$pipe1Name")
    updatePipeEndPath "$pipe1Name" "$tag"
    local r1_outpipe=$(getPipeEndPath "$pipe1Name")

    local r2_inpipe=$(getPipeEndPath "$pipe2Name")
    updatePipeEndPath "$pipe2Name" "$tag"
    local r2_outpipe=$(getPipeEndPath "$pipe2Name")

    # Bash SUCKS! Empty arrays do give an error with set -u with Bash < 4.4, but -v varName still succeeds! Therefore test the content, here.
    if [[ "${rest:-}" != "" ]]; then
        "$command" "$r1_inpipe" "$r2_inpipe" "$r1_outpipe" "$r2_outpipe" "${rest[@]}" & registerPid
    else
        "$command" "$r1_inpipe" "$r2_inpipe" "$r1_outpipe" "$r2_outpipe" & registerPid
    fi
}



#############################################################################

getFastqAsciiStream() {
    local name="${1:?No FASTQ filename given}"
    $UNZIPTOOL $UNZIPTOOL_OPTIONS "$name"
}

reorderUndirectionalReads() {
    local r1Input="${1:?No input R1}"
    local r2Input="${2:?No input R2}"
    local r1Output="${3:?No output R1}"
    local r2Output="${4:?No output R2}"
    "$PYTHON_BINARY" "$TOOL_UNIDIRECTIONAL_WGBS_READ_REORDERING" \
        --input_ascii \
        --R1_in "$r1Input" \
        --R2_in "$r2Input" \
        --output_ascii \
        --R1_out "$r1Output" \
        --R2_out "$r2Output" \
        --R1_unassigned /dev/null \
        --R2_unassigned /dev/null
}

toIlluminaScore() {
    local inFile="${1:?No input file}"
    local outFile="${2:?No output file}"
    "$PERL_BINARY" "$TOOL_CONVERT_ILLUMINA_SCORES" "$inFile" > "$outFile"
}

methylCfqconv() {
    local inFile="${1:?No input file}"
    local outFile="${2:?No output file}"
    local readNo="${3:?No read number}"
    "$PYTHON_BINARY" "$TOOL_METHYL_C_TOOLS" fqconv "-$readNo" "$inFile" "$outFile"
}

trimmomatic() {
    local input1="${1:?No R1 input}"
    local input2="${2:?No R2 input}"
    local output1="${3:?No R1 output}"
    local output2="${4:?No R2 output}"

    local u1=/dev/null
    local u2=/dev/null

    "$JAVA_BINARY" -jar "$TOOL_ADAPTOR_TRIMMING" \
        $ADAPTOR_TRIMMING_OPTIONS_0 \
        "$input1" "$input2" "$output1" "$u1" "$output2" "$u2" \
        $ADAPTOR_TRIMMING_OPTIONS_1
}


trimmomaticSingleEnd() {
    local input1="${1:?No R1 input}"
    local input2="${2:?No R2 input}"
    local output1="${3:?No R1 output}"
    local output2="${4:?No R2 output}"

    throw 254 "Not implemented"
}


preprocessWgbsPairedEndReads() {
    local fqName="${1:?No pipe-name given}"

    # Determine the quality score
    set +e   # The following few lines fail with exit 141 due to unknown reasons if `set -e` is set.
    qualityScore=$(getFastqAsciiStream "$RAW_SEQ_1" | "$PERL_BINARY" "$TOOL_SEQUENCER_DETECTION")
    set -e

    getFastqAsciiStream "$RAW_SEQ_1" > $(getPairedPipeEndPath 1 "$fqName") & registerPid
    getFastqAsciiStream "$RAW_SEQ_2" > $(getPairedPipeEndPath 2 "$fqName") & registerPid

    if [[ "$qualityScore" == "illumina" ]]; then
        extendPipe $(mkPairedPipeName 1 "$fqName") "qScore" -- toIlluminaScore
        extendPipe $(mkPairedPipeName 2 "$fqName") "qScore" -- toIlluminaScore
    fi

    if [[ "${useAdaptorTrimming:-false}" == "true" ]]; then
        extendPipePair "$fqName" "trimmomatic" -- trimmomatic
    fi

    if [[ "${reorderUndirectionalWGBSReadPairs:-false}" == "true" ]]; then
        if [[ "$useSingleEndProcessing" == "true" ]]; then
            throw 50 "Cannot reorder undirectional read pairs with single-end data"
        fi
        extendPipePair "$fqName" "reorder" -- reorderUndirectionalReads
    fi

    extendPipe $(mkPairedPipeName 1 "$fqName") "fqconv" -- methylCfqconv 1
	extendPipe $(mkPairedPipeName 2 "$fqName") "fqconv" -- methylCfqconv 2
}


preprocessWgbsSingleEndReads() {
    local fqName="${1:?No pipe-name given}"

    # Determine the quality score
    set +e   # The following few lines fail with exit 141 due to unknown reasons if `set -e` is set.
    qualityScore=$(getFastqAsciiStream "$RAW_SEQ_1" | "$PERL_BINARY" "$TOOL_SEQUENCER_DETECTION")
    set -e

    cat "$(getPairedPipeEndPath 2 $fqName)" > /dev/null
    getFastqAsciiStream "$RAW_SEQ_1" > $(getPairedPipeEndPath 1 "$fqName") & registerPid

    if [[ "$qualityScore" == "illumina" ]]; then
        extendPipe $(mkPairedPipeName 1 "$fqName") "qScore" -- toIlluminaScore
    fi

    if [[ "${useAdaptorTrimming:-false}" == "true" ]]; then
        extendPipePair "$fqName" "trimmomatic" -- trimmomaticSingleEnd
    fi

    if [[ "${reorderUndirectionalWGBSReadPairs:-false}" == "true" ]]; then
        if [[ "$useSingleEndProcessing" == "true" ]]; then
            throw 50 "Cannot reorder undirectional read pairs with single-end data"
        fi
    fi

    extendPipe $(mkPairedPipeName 1 "$fqName") "fqconv" -- methylCfqconv 1
}



eval "$WORKFLOWLIB___SHELL_OPTIONS"
if [[ "$WORKFLOWLIB___ERREXIT" == "errexit" ]]; then
    set -e
fi
