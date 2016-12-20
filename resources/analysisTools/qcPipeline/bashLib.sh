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
    hostname -f
    ulimit -a
    echo "user="$(whoami)
    echo "umask="$(umask)
    echo "groups="$(groups)
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
    declare -la values=($@)
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


eval "$BASHLIB___SHELL_OPTIONS"

