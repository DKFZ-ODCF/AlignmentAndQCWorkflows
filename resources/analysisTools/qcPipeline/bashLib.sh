# Library of BASH function. Please import using
#
#        source "$TOOL_BASH_LIB"


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
    {
        hostname -f
        ulimit -a
        echo "user="$(whoami)
        echo "umask="$(umask)
        echo "groups="$(groups)
    } >> /dev/stderr
}


toMinusIEqualsList () {
    declare -a inputFiles=($@)
    for inFile in ${inputFiles[@]}; do
            echo -n "I=$inFile "
    done
    echo
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