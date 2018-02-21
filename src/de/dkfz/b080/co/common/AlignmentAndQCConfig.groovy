package de.dkfz.b080.co.common

import de.dkfz.b080.co.files.COConstants
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

    public AlignmentAndQCConfig(ExecutionContext context) {
        super(context)
    }

    public String getSingleBamParameter() {
        return configValues.get("bam", "");
    }

    boolean getUseOnlyExistingTargetBam() {
        return configValues.getBoolean("useOnlyExistingTargetBam", false)
    }

    boolean getUseExistingLaneBams() {
        return configValues.getBoolean(COConstants.FLAG_USE_EXISTING_LANE_BAMS, false)
    }

    public String getIndexPrefix() {
        return configValues.getString(CVALUE_INDEX_PREFIX, "")
    }

    public File getChromosomeSizesFile() {
        return new File (configValues.getString(CVALUE_CHROMOSOME_SIZES_FILE, ""))
    }

    public File getTargetRegionsFile() {
        return new File (configValues.getString(CVALUE_TARGET_REGIONS_FILE, ""))
    }

    public Integer getTargetSize() {
        Integer returnValue = configValues.getString(CVALUE_TARGET_SIZE, null) as Integer
        if (null == returnValue) {
            returnValue = configValues.getString(CVALUE_TARGETSIZE, null) as Integer
        }
        return returnValue
    }

    public boolean getRunExomeAnalysis() {
        return configValues.getBoolean(COConstants.FLAG_RUN_EXOME_ANALYSIS)
    }

    public File getCytosinePositionIndex() {
        return new File(configValues.getString(CVALUE_CYTOSINE_POSITIONS_INDEX))
    }

    public File getClipIndex() {
        return new File(configValues.getString(CVALUE_CLIP_INDEX))
    }

    public Boolean getUseOnlyExistingPairedBams() {
        if (configValues.hasValue(COConstants.FLAG_USE_EXISTING_LANE_BAMS)) {
            System.err.println("${COConstants.FLAG_USE_EXISTING_LANE_BAMS} is deprecated. Use ${COConstants.FLAG_USE_ONLY_EXISTING_PAIRED_BAMS}.")
            return configValues.getBoolean(COConstants.FLAG_USE_EXISTING_LANE_BAMS, false);
        } else {
            return configValues.getBoolean(COConstants.FLAG_USE_ONLY_EXISTING_PAIRED_BAMS, false);
        }
    }

}
