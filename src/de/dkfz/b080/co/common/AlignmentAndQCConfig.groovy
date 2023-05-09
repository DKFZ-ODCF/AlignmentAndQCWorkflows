package de.dkfz.b080.co.common

import de.dkfz.b080.co.files.COConstants
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.core.ExecutionContext
import groovy.transform.CompileStatic

/**
 * Created by kensche on 17.05.16.
 */
@CompileStatic
class AlignmentAndQCConfig extends COConfig {

    public static final String CVALUE_INDEX_PREFIX = "INDEX_PREFIX"
    public static final String CVALUE_CHROMOSOME_SIZES_FILE = "CHROM_SIZES_FILE"
    public static final String CVALUE_TARGET_REGIONS_FILE = "TARGET_REGIONS_FILE"
    public static final String CVALUE_TARGETSIZE = "TARGETSIZE"
    public static final String CVALUE_TARGET_SIZE = "TARGET_SIZE"
    public static final String CVALUE_CLIP_INDEX = "CLIP_INDEX"
    public static final String CVALUE_CYTOSINE_POSITIONS_INDEX = "CYTOSINE_POSITIONS_INDEX"
    public static final String CVALUE_RUN_FINGERPRINTING = "runFingerprinting"
    public static final String CVALUE_FINGERPRINTING_SITES_FILE="fingerprintingSitesFile"

    AlignmentAndQCConfig(ExecutionContext context) {
        super(context)
    }

    void setUseOnlyExistingTargetBam(boolean value = true) {
        setConfig("useOnlyExistingTargetBam", value.toString(), "boolean")
    }

    void setExtractSamplesFromOutputFiles(boolean value = true) {
        setConfig("extractSamplesFromOutputFiles", value.toString(), "boolean")
    }

    String getSingleBamParameter() {
        return configValues.get("bam", "");
    }

    boolean getUseOnlyExistingTargetBam() {
        return configValues.getBoolean("useOnlyExistingTargetBam", false)
    }

    boolean getUseExistingLaneBams() {
        return configValues.getBoolean("useExistingLaneBams", false)
    }

    boolean getUseOnlyExistingPairedBams() {
        return configValues.getBoolean("useExistingPairedBams", false)
    }

    String getIndexPrefix() {
        return configValues.getString(CVALUE_INDEX_PREFIX, "")
    }

    File getChromosomeSizesFile() {
        return new File (configValues.getString(CVALUE_CHROMOSOME_SIZES_FILE, ""))
    }

    File getTargetRegionsFile() {
        return new File (configValues.getString(CVALUE_TARGET_REGIONS_FILE, ""))
    }

    Integer getTargetSize() {
        Integer returnValue = configValues.getString(CVALUE_TARGET_SIZE, null) as Integer
        if (null == returnValue) {
            returnValue = configValues.getString(CVALUE_TARGETSIZE, null) as Integer
        }
        return returnValue
    }

    boolean getRunExomeAnalysis() {
        return configValues.getBoolean("runExomeAnalysis")
    }

    File getCytosinePositionIndex() {
        return new File(configValues.getString(CVALUE_CYTOSINE_POSITIONS_INDEX))
    }

    File getClipIndex() {
        return new File(configValues.getString(CVALUE_CLIP_INDEX))
    }

    boolean getRunFingerprinting() {
        return configValues.getBoolean(CVALUE_RUN_FINGERPRINTING, true)
    }

    File getFingerprintingSitesFile() {
        return new File(configValues.getString(CVALUE_FINGERPRINTING_SITES_FILE))
    }

    boolean getRunFastqcOnly() {
        return configValues.getBoolean("runFastQCOnly", false)
    }

    boolean getRunFastqc() {
        return configValues.getBoolean("runFastQC", true)
    }

    boolean getRunAlignmentOnly() {
        return configValues.getBoolean("runAlignmentOnly", false)
    }

    boolean getRunCoveragePlots() {
        return configValues.getBoolean("runCoveragePlots", true)
    }

    boolean getRunCollectBamFileMetrics() {
        return configValues.getBoolean("runCollectBamFileMetrics", false)
    }

    List<String> getOverrideSampleNames() {
        return configValues.getString("overrideSampleNames", "").
                split(StringConstants.SPLIT_SEMICOLON) as List<String>
    }

    @Deprecated
    boolean getUseCombinedAlignAndSampe() {
        return true
    }

    @Deprecated
    boolean getRunSlimWorkflow() {
        return true
    }

}
