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
    assertFalse "arrayContains x"
    assertFalse "arrayContains x a b c d"
    assertTrue "arrayContains x x b c d"
    assertTrue "arrayContains x a b x d"
    assertTrue "arrayContains x a b c x"
}

testMatchPrefixInArray() {
    assertFalse "matchPrefixInArray x"
    assertFalse "matchPrefixInArray x a b cx"
    assertTrue "matchPrefixInArray x x a b c"
    assertTrue "matchPrefixInArray x xa a b c"
    assertTrue "matchPrefixInArray x a xb c"
    assertTrue "matchPrefixInArray x a b xc"
}

testIsControlSample() {
    local possibleControlSampleNamePrefixes="(ax by cz)"
    assertFalse "isControlSample x"
    assertFalse "isControlSample z"
    assertTrue "isControlSample a"
    assertTrue "isControlSample b"
    assertTrue "isControlSample c"
    assertTrue "isControlSample cz"
}

testIsTumorSample() {
    local possibleTumorSampleNamePrefixes="(ax by cz)"
    assertFalse "isTumorSample x"
    assertFalse "isTumorSample z"
    assertTrue "isTumorSample a"
    assertTrue "isTumorSample b"
    assertTrue "isTumorSample c"
    assertTrue "isTumorSample cz"
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

testMatchesShortChromosomeName() {
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
    local CHROMOSOME_INDICES="(1 2 3 chrMmu1 chrMmuX)"
    local CHR_PREFIX=chrMmu

    assertEquals "short=1,2,3" $(shortChromosomeGroupSpec)

    SHORT_CHROMOSOME_NAME_GROUP=human
    assertEquals "human=1,2,3" $(shortChromosomeGroupSpec)

}

testLongChromosomeGroupSpec() {
    local CHROMOSOME_INDICES="(1 2 3 chrMmu1 chrMmuX)"
    local CHR_PREFIX=chrMmu

    assertEquals "long=chrMmu1,chrMmuX" $(longChromosomeGroupSpec)

    LONG_CHROMOSOME_NAME_GROUP=mouse
    assertEquals "mouse=chrMmu1,chrMmuX" $(longChromosomeGroupSpec)
}

source ${SHUNIT2:?Oops}