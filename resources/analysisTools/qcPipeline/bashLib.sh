# Library of BASH function. Please import using
#
#        source "$TOOL_BASH_LIB"


BASHLIB___SHELL_OPTIONS=$(set +o)
set +o verbose
set +o xtrace

UNSPECIFIED_ERROR_MSG="Unspecified error"
UNSPECIFIED_ERROR_CODE=1

EMPTY_VALUE_MSG="Empty value"
EMPTY_VALUE_CODE=200

NOT_WRITABLE_MSG="File or directory is not writable"
NOT_WRITABLE_CODE=201

## From http://unix.stackexchange.com/questions/26676/how-to-check-if-a-shell-is-login-interactive-batch
shellIsInteractive () {
    case $- in
        *i*) echo "true";;
        *)   echo "false";;
    esac
}


## funname () ( set +exv; ...; ) may be better to get rid of too much output (mind the (...) subshell) but the exit won't work anymore.
## Maybe set -E + trap "bla" ERR would work? http://fvue.nl/wiki/Bash:_Error_handling#Exit_on_error
printStackTrace () {
    frameNumber=0
    while caller $frameNumber ;do
      ((frameNumber++))
    done
}


errout () {
    local exitCode="$1"
    local message="$2"
    env printf "Error(%d): %s\n" "$exitCode" "$message" >> /dev/stderr
}


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


printInfo () {
    ## Get information about the node.
    {
        hostname -f
        ulimit -a
        echo "user="$(whoami)
        echo "umask="$(umask)
        echo "groups="$(groups)
    } >> /dev/stderr
}


## Use 'assertNonEmpty $parameter || return $?'
assertNonEmpty () {

    local value="$1"
    local message="${2-EMPTY_VALUE_MSG}"
    if [[ "$value" == "" ]]; then
        throw "$EMPTY_VALUE_CODE" "$message" || return $?
    fi
}


waitAndMaybeExit () {
    local pid="$1"
    local errorMessage="${2-$UNSPECIFIED_ERROR_MSG}"
    local errorCode="${3-$UNSPECIFIED_ERROR_CODE}"
    wait $pid; returnValue=$?
    if [[ $returnValue -gt 0 ]]; then
        throw "$errorCode" "$errorMessage" || return $?
    fi
}


stringJoin () {
    local separator="$1"
    shift
    assertNonEmpty "$separator" "Undefined separator" || return $?
    declare -a values=($@)
    local result=""
    local first=true
    for value in ${values[@]}; do
        if [[ $first == true ]]; then
            result="$value"
            first=false
        else
            result="$result$separator$value"
        fi
    done
    echo "$result"
}

normalizeBoolean() {
    if [[ "${1:-false}" == "true" ]]; then
        echo "true"
    else
        echo "false"
    fi
}

isDebugSet() {
    normalizeBoolean "${debug:-false}"
}


#####################################################################
## Handling processes and tempfiles (original code from BamToFastqPlugin)
#####################################################################

# TODO: Make this based on associative array to have more information (some tag!) about the failed process. Every tag only once (warn and save with extension!)
registerPid() {
    local pid="${1:-$!}"
    declare -gax pids=(${pids[@]} $pid)
}

listPids() {
    for pid in "${pids[@]}"; do
        if [[ "$pid" != "$ARRAY_ELEMENT_DUMMY" ]]; then
            echo "$pid"
        fi
    done
}

registerTmpFile() {
    local tmpFile="${1:?No temporary file name to register}"
    # Note that the array is build in reversed order, which simplifies the deletion of nested directories.
    declare -gax tmpFiles=("$tmpFile" "${tmpFiles[@]}")
}

reverseArray() {
    local c=""
    for b in "$@"; do
        c="$b $c"
    done
    echo $c
}

# Bash sucks. An empty array does not exist! So if there are no tempfiles/pids, then there is no array and set -u will result an error!
# Therefore we put a dummy value into the arrays and have to take care to remove the dummy before the processing.
# The dummy contains a random string to avoid collision with possible filenames (the filename 'dummy' is quite likely).
ARRAY_ELEMENT_DUMMY=$(mktemp -u "_dummy_XXXXX")

waitForRegisteredPids_BashSucksVersion() {
    jobs
    declare -a realPids=$(listPids)
    if [[ -v realPids && ${#realPids[@]} -gt 0 ]]; then
        wait ${realPids[@]}
        declare EXIT_CODE=$?
        if [[ ${EXIT_CODE} -ne 0 ]]; then
            throw ${EXIT_CODE} "One of the following processes ended with exit code ${EXIT_CODE}: ${realPids[@]}"
        fi
    fi
    pids=("$ARRAY_ELEMENT_DUMMY")
}
setUp_BashSucksVersion() {
    declare -g -a -x tmpFiles=("$ARRAY_ELEMENT_DUMMY")
    declare -g -a -x pids=("$ARRAY_ELEMENT_DUMMY")
    initPipeEnds

    # Remove all registered temporary files upon exit
    trap cleanUp_BashSucksVersion EXIT
}
cleanUp_BashSucksVersion() {
    if [[ $(isDebugSet) == "false" && -v tmpFiles && ${#tmpFiles[@]} -gt 1 ]]; then
        for f in ${tmpFiles[@]}; do
            if [[ "$f" == "$ARRAY_ELEMENT_DUMMY" ]]; then
                continue
            elif [[ -d "$f" ]]; then
                rmdir "$f"
            elif [[ -e "$f" ]]; then
                rm -f "$f"
            fi
        done
        tmpFiles=("$ARRAY_ELEMENT_DUMMY")
    fi
}

# These versions only works with Bash >4.4. Prior version do not really declare the array variables with empty values and set -u results in error message.
waitForRegisteredPids() {
    jobs
    wait ${pids[@]}
    pids=()
}
setUp() {
    declare -g -a -x tmpFiles=()
    declare -g -a -x pids=()
    initPipeEnds
}
cleanUp() {
    if [[ $(isDebugSet) == "false" && -v tmpFiles && ${#tmpFiles[@]} -gt 0 ]]; then
        for f in "${tmpFiles[@]}"; do
            if [[ -d "$f" ]]; then
                rmdir "$f"
            elif [[ -e "$f" ]]; then
                rm "$f"
            fi
        done
        tmpFiles=()
    fi
}



#####################################################################
## Linear pipe management
#####################################################################
# The pipe-extension API follows this pattern:
#
# * pipeExtenderFunction pipeBaseName pipeStepExtension -- command
#
# * The pipeBaseName is the name under which the pipe (or derived from this: pipe pair) is registered.
# * The pipeStepExtension is used to tag the actual pipe used in the specific step. Use a name describing the content flowing through the pipe,
#   i.e. the output of the extension step you are declaring.
# * Current pipe extender functions are extendPipe (for single pipes) and extendPipePair (for r1/r2 pipes)
# * To work with paired pipes (r1/r2) mkPairedPipeName is used to calculate the r1_ or r2_ pipeBaseName

_pipePath="$RODDY_SCRATCH"

# Maintain a mapping of pipe-basenames to current pipe-ends.
initPipeEnds() {
    unset _pipeEnds
    declare -Ag _pipeEnds=()
}

listPipeEnds() {
    for i in "${!_pipeEnds[@]}"; do
        echo "'$i' => '${_pipeEnds[$i]}'"
    done
}

# Make a path without registering it.
mkPipePath() {
    local pipeName="${1:?Named-pipe base name}"
    echo "$_pipePath/$pipeName"
}

# Just get the path from the registry.
getPipeEndPath() {
    local pipeName="${1:?No pipename}"
    if [[ ! ${_pipeEnds[$pipeName]+_} ]]; then
        throw 152 "Could not find pipe by name '$pipeName'"
    fi
    local pipePath="${_pipeEnds[$pipeName]}"
    echo "$pipePath"
}

# Register a path in the registry under the key "pipeName".
setPipeEndPath() {
    local pipeName="${1:?No pipename}"
    local pipePath="${2:?No pipe path}"
    if [[ -z "$pipeName" ]]; then
        throw 153 "Cannot set empty pipename"
    fi
    _pipeEnds[$pipeName]="$pipePath"
}

# Create a path from the name and the tag and register it as new pipe-end in the registry.
# Note that the pipe is registered as tempfile and will be deleted upon cleanUp.
updatePipeEndPath() {
    local pipeName="${1:?Named-pipe base name}"
    local tag="${2:?Pipe tag}"
    local pipePath=$(mkPipePath "${pipeName}_$tag")
    mkfifo "$pipePath" || throw 150 "Could not create named pipe at '$pipePath'"
    registerTmpFile "$pipePath"
    setPipeEndPath "$pipeName" "$pipePath"
}

# Create a pipe and register it as new pipe-end in the registry. It will be the source of pipe.
mkPipeSource() {
    local pipeName="${1:?No pipename}"
    if [[ ${_pipeEnds[$pipeName]+_} ]]; then
        throw 154 "Cannot create new pipe source. Pipe named '$pipeName' -- already exists with value  '${_pipeEnds[$pipeName]}'"
    fi
    updatePipeEndPath "$pipeName" "source"
}

# Linear extension of a pipeline.
# For a command with the interface "command infile outfile @rest" take the pipe with the given basename
# and use it as input file. Create a new output pipe and set that pipe as output file.
# Set the new end of the linear pipeline to the output pipe. Call like this:
#
#   extendPipe $pipeName $tag -- commandName @restArgs
#
# The "--" is optional.
extendPipe() {
    local pipeName="${1:?No pipe basename}"
    local tag="${2:?No tag}"
    shift 2
    if [[ "$1" == "--" ]]; then
        shift
    fi
    local command="${1:?No command/function}"
    shift
    declare -a args=("$@")

    local inpipe=$(getPipeEndPath "$pipeName")
    updatePipeEndPath "$pipeName" "$tag"
    local outpipe=$(getPipeEndPath "$pipeName")

    "$command" "$inpipe" "$outpipe" "${args[@]}" & registerPid
}


eval "$BASHLIB___SHELL_OPTIONS"

