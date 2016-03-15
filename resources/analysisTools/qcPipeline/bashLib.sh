# Library of BASH function. Please import using
#
#        source "$TOOL_BASH_LIB"


ILLEGAL_ARGUMENT_ERROR=200

## throw [code [msg]]
## Write message (Unspecified error) to STDERR and exit with code (default 1)
function throw {
  local exitCode="${1-1}"
  local msg="${2-Unspecified error}"
  echo "Error($exitCode): file $BASH_SOURCE, function ${FUNCNAME[1]}, line ${BASH_LINENO[0]}: $msg" >> /dev/stderr
  exit "$exitCode"
}

function printInfo {
    ## Get information about the node.
    hostname -f
    ulimit -a
    echo "user="$(whoami)
    echo "umask="$(umask)
    echo "groups="$(groups)
}

function runningOnConvey {
    if [[ "$PBS_QUEUE" == convey* ]]; then
    	echo "true"
    else
    	echo "false"
    fi
}

# Stolen from http://stackoverflow.com/questions/3685970/check-if-an-array-contains-a-value
function arrayContains {
    local ELEMENT="${1}"
    if [[ "$ELEMENT" == "" ]]; then
        throw "arrayContains called without parameters"
    fi
    local DELIM=","
    printf "${DELIM}%s${DELIM}" "${@:2}" | grep -q "${DELIM}${ELEMENT}${DELIM}"
}

function matchPrefixInArray {
    local ELEMENT="${1}"
    if [[ "$ELEMENT" == "" ]]; then
        throw "matchPrefixInArray called without parameters"
    fi
    local DELIM=","
    printf "${DELIM}%s${DELIM}" "${@:2}" | grep -q -P "${DELIM}${ELEMENT}[^${DELIM}]*${DELIM}"
}

function isControlSample {
    if [[ "$1" == "" ]]; then
        throw $ILLEGAL_ARGUMENT_ERROR "isControlSample expects sample type name as single parameter"
    fi
    if [[ "$possibleControlSampleNamePrefixes" == "" ]]; then
        throw $ILLEGAL_ARGUMENT_ERROR "Undefined/empty possibleControlSampleNamePrefixes"
    fi
    matchPrefixInArray "$1" "${possibleControlSampleNamePrefixes[@]}"
}

function isTumorSample {
    if [[ "$1" == "" ]]; then
        throw $ILLEGAL_ARGUMENT_ERROR "isTumorSample expects sample type name as single parameter"
    fi
    if [[ "$possibleTumorSampleNamePrefixes" == "" ]]; then
        throw $ILLEGAL_ARGUMENT_ERROR "Undefined/empty possibleTumorSampleNamePrefixes"
    fi
    matchPrefixInArray "$1" "${possibleTumorSampleNamePrefixes[@]}"
}

function sampleType {
    if [[ "$1" == "" ]]; then
        throw $ILLEGAL_ARGUMENT_ERROR "sampleType expects sample type name as single parameter"
    fi
    if isControlSample "$1"; then
        echo "control"
    elif isTumorSample "$1"; then
        echo "tumor"
    else
        throw $ILLEGAL_ARGUMENT_ERROR "'$1' is neither control nor tumor"
    fi
}
