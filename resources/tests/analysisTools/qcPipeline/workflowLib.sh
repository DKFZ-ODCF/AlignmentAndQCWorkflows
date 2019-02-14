#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

source ${TOOL_WORKFLOW_LIB:?No TOOL_WORKFLOW_LIB}


testMarkWithPicard() {
    assertFalse 1 markWithPicard

    local markDuplicatesVariant=picard
    assertTrue 2 markWithPicard
    local markDuplicatesVariant=sambamba
    assertFalse 3 markWithPicard

    unset markDuplicatesVariant
    local useBioBamBamMarkDuplicates=false
    assertTrue 4 markWithPicard
}

testMarkWithSambamba() {
    assertFalse markWithSambamba
    local markDuplicatesVariant=sambamba
    assertTrue markWithSambamba
}

testMarkWithBiobambam() {
    assertFalse markWithBiobambam

    local markDuplicatesVariant=biobambam
    assertTrue markWithBiobambam
    local useBiobambamMarkDuplicates=false
    assertTrue markWithBiobambam
}


testGetBigScratchDirectory() {
    local outputAnalysisBaseDirectory=/tmp/blub
    local RODDY_SCRATCH=/tmp/blabla

    assertEquals $outputAnalysisBaseDirectory/tmp "$(getBigScratchDirectory)"
    assertEquals /tmp/alt "$(getBigScratchDirectory /tmp/alt)"

    local useRoddyScratchAsBigFileScratch=true
    assertEquals $RODDY_SCRATCH "$(getBigScratchDirectory)"
}

testAnalysisType() {
    assertEquals genome "$(analysisType)"
    local runExomeAnalysis=true
    assertEquals exome "$(analysisType)"
}


testToIEqualsList() {
    assertEquals "I=a I=b I=c " "$(toIEqualsList a b c)"
}


setupPipePath() {
    initPipeEnds
    _pipePath=$(mktemp -d /tmp/bashlibTest_XXXXXX)
}

teardownPipePath() {
    rm -rf "$_pipePath"
    unset _pipePath
    unset _pipeEnds
}


testMkPairedPipeName() {
    local name

    name=$(mkPairedPipeName 1 tag)
    assertEquals "r1_tag" "$name"

    name=$(mkPairedPipeName 2 tag)
    assertEquals "r2_tag" "$name"
}

testMkPipePair() {
    setupPipePath
    mkPipePairSource "test"
    test -p $(getPairedPipeEndPath 1 "test")
    assertTrue "First pipe in pair created" $?
    test -p $(getPairedPipeEndPath 2 "test")
    assertTrue "Second pipe in pair created" $?
    teardownPipePath
}

reorder() {
    cat "$1" <(echo "_was1") > "$4" &
    cat "$2" <(echo "_was2") > "$3" &
    wait
}

testExtendPipePair() {
    setupPipePath
    mkPipePairSource "test"

    local source1=$(getPairedPipeEndPath 1 "test")
    echo -n "hallo1" > "$source1" &

    local source2=$(getPairedPipeEndPath 2 "test")
    echo -n "hallo2" > "$source2" &

    extendPipePair "test" "step1" -- reorder

    local result1=$(cat $(getPairedPipeEndPath 1 "test"))
    local result2=$(cat $(getPairedPipeEndPath 2 "test"))

    assertEquals "hallo1_was1" "$result2"
    assertEquals "hallo2_was2" "$result1"

    teardownPipePath
}

source ${SHUNIT2:?Oops}