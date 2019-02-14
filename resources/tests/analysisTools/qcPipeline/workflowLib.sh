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

testArrayContains() {
    assertFalse "Fail on empty arrai" "arrayContains x"
    assertFalse "Fail non-existing" "arrayContains x a b c d"
    assertTrue "Succeed existing at start" "arrayContains x x b c d"
    assertTrue "Succeed existing in the middle" "arrayContains x a b x d"
    assertTrue "Succeed existing in the end" "arrayContains x a b c x"
}

testMatchPrefixInArray() {
    assertFalse "Empty match prefix" "matchPrefixInArray x"
    assertFalse "Empty match prefix with suffix" "matchPrefixInArray x a b cx"
    assertTrue "Match fully" "matchPrefixInArray x x a b c"
    assertTrue "Match prefix at start" "matchPrefixInArray x xa a b c"
    assertTrue "Match prefix in the middle" "matchPrefixInArray x a xb c"
    assertTrue "Match prefix in the end" "matchPrefixInArray x a b xc"
}

testIsControlSample() {
    local possibleControlSampleNamePrefixes="(ax by cz)"

    assertFalse "Mismatch control sample 1" "isControlSample x"
    assertFalse "Mismatch control sample 2" "isControlSample z"
    assertTrue "Match control prefix beginning" "isControlSample a"
    assertTrue "Match control prefix middle" "isControlSample b"
    assertTrue "Match control prefix end" "isControlSample c"
    assertTrue "Match control fully" "isControlSample cz"
}

testIsTumorSample() {
    local possibleTumorSampleNamePrefixes="(ax by cz)"

    assertFalse "Mismatch tumor 1" "isTumorSample x"
    assertFalse "Mismatch tumor 2" "isTumorSample z"
    assertTrue "Match control prefix start" "isTumorSample a"
    assertTrue "Match control prefix middle" "isTumorSample b"
    assertTrue "Match control prefix end" "isTumorSample c"
    assertTrue "Match control fully" "isTumorSample cz"
}

testSampleType() {
    declare -x possibleControlSampleNamePrefixes="(cA cB vZ)"
    declare -x possibleTumorSampleNamePrefixes="(tA tB vZ)"

    declare -fx assertNonEmpty sampleType isControlSample isTumorSample throw_illegal_argument throw matchPrefixInArray errout printStackTrace exitIfNonInteractive shellIsInteractive
    assertFalse "bash -c "sampleType v" 2> /dev/null"
    assertFalse "bash -c "sampleType vZ" 2> /dev/null"
    assertFalse "bash -c "sampleType q" 2> /dev/null"
    assertEquals control $(bash -c 'sampleType c')
    assertEquals control $(bash -c 'sampleType c')
    assertEquals tumor $(bash -c 'sampleType t')
    assertEquals tumor $(bash -c 'sampleType tA')
}

chromosomeSizesFile() {
    echo "$(dirname ${BASH_SOURCE[0]})/chrom-sizes-file.tsv"
}

testChromosomeIndices() {
    local CHROM_SIZES_FILE=$(chromosomeSizesFile)
    assertEquals "1 2 3 chrMmu1 chrMmuX" "$(chromosomeIndices)"
}

testMatchesShortChromosomeName() {
    assertEquals "true" $(matchesShortChromosomeName 1)
    assertEquals "true" $(matchesShortChromosomeName chr1)

    local CHR_PREFIX=chrMmu
    assertEquals "true"  $(matchesShortChromosomeName 1)
    assertEquals "true"  $(matchesShortChromosomeName 1xxx)
    assertEquals "false" $(matchesShortChromosomeName chrMmu1)
    assertEquals "false" $(matchesShortChromosomeName chrMmu1xxx)

    local CHR_SUFFIX=bla
    assertEquals "true"  $(matchesShortChromosomeName 1)
    assertEquals "true"  $(matchesShortChromosomeName 1xxx)
    assertEquals "true"  $(matchesShortChromosomeName 1bla)
    assertEquals "true"  $(matchesShortChromosomeName chrMmu1)
    assertEquals "true"  $(matchesShortChromosomeName chrMmu1xxx)
    assertEquals "false" $(matchesShortChromosomeName chrMmu1bla)
}


testMatchesLongChromosomeName() {
    assertEquals "false" $(matchesLongChromosomeName 1)
    assertEquals "false" $(matchesLongChromosomeName chr1)

    local CHR_PREFIX=chrMmu
    assertEquals "false" $(matchesLongChromosomeName 1)
    assertEquals "false" $(matchesLongChromosomeName 1xxx)
    assertEquals "true"  $(matchesLongChromosomeName chrMmu1)
    assertEquals "true"  $(matchesLongChromosomeName chrMmu1xxx)

    local CHR_SUFFIX=bla
    assertEquals "false" $(matchesLongChromosomeName 1)
    assertEquals "false" $(matchesLongChromosomeName 1xxx)
    assertEquals "false" $(matchesLongChromosomeName 1bla)
    assertEquals "false" $(matchesLongChromosomeName chrMmu1)
    assertEquals "false" $(matchesLongChromosomeName chrMmu1xxx)
    assertEquals "true"  $(matchesLongChromosomeName chrMmu1bla)
}

testShortChromosomeGroupSpec() {
    local CHROM_SIZES_FILE=$(chromosomeSizesFile)
    local CHR_PREFIX=chrMmu

    assertEquals "nonmatching=1,2,3" $(shortChromosomeGroupSpec)

    CHR_GROUP_NOT_MATCHING=human
    assertEquals "human=1,2,3" $(shortChromosomeGroupSpec)

}

testLongChromosomeGroupSpec() {
    local CHROM_SIZES_FILE=$(chromosomeSizesFile)
    local CHR_PREFIX=chrMmu

    assertEquals "matching=chrMmu1,chrMmuX" $(longChromosomeGroupSpec)

    CHR_GROUP_MATCHING=mouse
    assertEquals "mouse=chrMmu1,chrMmuX" $(longChromosomeGroupSpec)
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