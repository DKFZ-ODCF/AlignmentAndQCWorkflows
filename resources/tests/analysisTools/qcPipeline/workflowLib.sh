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

source ${SHUNIT2:?Oops}