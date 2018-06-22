#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
#

source ${TOOL_BASH_LIB:?No TOOL_BASHLIB}

testShellIsInteractive() {
    export -f shellIsInteractive
    assertEquals false "$(bash -c shellIsInteractive)"
    assertEquals true "$(bash -i -c shellIsInteractive)"
}

testPrintErrout() {
    assertEquals "Error(1): hola" "$(2>&1 errout 1 'hola')"
}

testAssertNonEmpty() {
    assertEquals "" "$(bash -i -c 'assertNonEmpty 2> /dev/null')"
    assertEquals "" "$(assertNonEmpty 'hola')"
}

testStringJoin() {
    assertEquals "" "$(stringJoin ',')"
    assertEquals "A:b:c" "$(stringJoin ':' A b c)"
    assertEquals "a,B" "$(stringJoin ',' a B)"
}

source ${SHUNIT2:?Oops}