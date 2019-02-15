##############################################################################
## Domain-specific code (maybe put this into a dedicated library file)
##############################################################################
##
## Include into your code with: source "$TOOL_WORKFLOW_LIB"
##

source "$TOOL_BASH_LIB"

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
    "$MBUFFER_BINARY" -m "$bufferSize" -q -l /dev/null ${@}
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
    local result
    result=$("$TOOL_BAM_IS_COMPLETE" "$bamFile")
    if [[ $? ]]; then
        echo "BAM is terminated! $bamFile" >> /dev/stderr
    else
        throw 40 "BAM is not terminated! $bamFile"
    fi
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
    local args=("$@")

    local pipe1Name=$(mkPairedPipeName 1 "$pipeName")
    local pipe2Name=$(mkPairedPipeName 2 "$pipeName")

    local r1_inpipe=$(getPipeEndPath "$pipe1Name")
    updatePipeEndPath "$pipe1Name" "$tag"
    local r1_outpipe=$(getPipeEndPath "$pipe1Name")

    local r2_inpipe=$(getPipeEndPath "$pipe2Name")
    updatePipeEndPath "$pipe2Name" "$tag"
    local r2_outpipe=$(getPipeEndPath "$pipe2Name")

    "$command" "$r1_inpipe" "$r2_inpipe" "$r1_outpipe" "$r2_outpipe" "${args[@]}" & registerPid
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
    "$TOOL_UNIDIRECTIONAL_WGBS_READ_REORDERING" \
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
    local i1="${1:?No R1 input}"
    local i2="${2:?No R2 input}"
    local o1="${3:?No R1 output}"
    local o2="${4:?No R2 output}"

    local u1=/dev/null
    local u2=/dev/null

    "$TRIMMOMATIC_BINARY" "$ADAPTOR_TRIMMING_OPTIONS_0" "$i1" "$i2" "$o1" "$u1" "$o2" "$u2" $ADAPTOR_TRIMMING_OPTIONS_1
}



eval "$WORKFLOWLIB___SHELL_OPTIONS"