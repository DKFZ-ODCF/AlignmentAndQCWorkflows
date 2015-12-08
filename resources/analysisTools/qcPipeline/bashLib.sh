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