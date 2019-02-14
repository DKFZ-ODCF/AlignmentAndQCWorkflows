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

source ${SHUNIT2:?Oops}