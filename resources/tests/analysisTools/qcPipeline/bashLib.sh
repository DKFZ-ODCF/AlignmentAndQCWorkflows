#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
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



testMkPipePath() {
    assertEquals "$_pipePath/test" "$(mkPipePath 'test')"
}

setupPipePath() {
    declare -A _pipeEnds=()
    _pipePath=$(mktemp -d /tmp/bashlibTest_XXXXXX)
}

teardownPipePath() {
    rm -rf "$_pipePath"
    unset _pipePath
    unset _pipeEnds
}

testSetGetPipeEndPath() {
    setupPipePath
    assertFalse "test getPipeEndPath 'doesnotexist'"
    setPipeEndPath "test" "/a/b/c"
    assertEquals "/a/b/c" "$(getPipeEndPath 'test')"
    teardownPipePath
}

testUpdatePipeEndPath() {
    setupPipePath
    updatePipeEndPath "test" "tag1"
    local sourcePath1=$(getPipeEndPath 'test')
    updatePipeEndPath "test" "tag2"
    local sourcePath2=$(getPipeEndPath 'test')
    assertNotEquals "$sourcePath1" "$sourcePath2"
    test -p "$sourcePath1"
    assertTrue "$sourcePath1 exists" $?
    test -p "$sourcePath2"
    assertTrue "$sourcePath2 exists" $?
    assertEquals "$sourcePath2" "$(getPipeEndPath 'test')"
    teardownPipePath
}

testMkPipeSource() {
    setupPipePath
    mkPipeSource "test"
    test -p $(getPipeEndPath 'test')
    assertTrue "source pipe end exists" $?
    #teardownPipePath
}

source ${SHUNIT2:?Oops}