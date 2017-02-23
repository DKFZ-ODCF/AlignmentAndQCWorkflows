package de.dkfz.b080.co.common

import de.dkfz.b080.co.files.COConstants
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextLevel
import de.dkfz.roddy.execution.io.ExecutionService
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

    public AlignmentAndQCConfig(ExecutionContext context) {
        super(context)
    }

    protected String substituteDirExecution(String filename) {
        return filename.replaceAll("[\$]${ExecutionService.RODDY_CVALUE_DIRECTORY_EXECUTION}(?=\\W|\$)", dirExecution)
    }

    public String getSingleBamParameter() {
        return configValues.get("bam", "");
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

    public boolean getUseAdapterTrimming() {
        return configValues.getBoolean(COConstants.FLAG_USE_ADAPTOR_TRIMMING, false)
    }

    public File getClipIndex() {
        return new File(configValues.getString(CVALUE_CLIP_INDEX))
    }

    public boolean getRunFingerprinting() {
        return configValues.getBoolean(CVALUE_RUN_FINGERPRINTING, true)
    }

    public File getFingerprintingSitesFile() {
        return new File (configValues.getString(CVALUE_FINGERPRINTING_SITES_FILE))
    }

    public String getDirExecution() {
        return context.executionDirectory.absolutePath
    }

    public boolean getRunFastQCOnly() {
        return configValues.getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false)
    }

    public boolean getRunAlignmentOnly() {
        return configValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false)
    }

    public boolean getRunCoveragePlots() {
        return configValues.getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true)
    }

    public boolean getRunSlimWorkflow() {
        return configValues.getBoolean(COConstants.FLAG_RUN_SLIM_WORKFLOW, false)
    }

    public boolean getRunCollectBamFileMetrics() {
        return configValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false)
    }

    public boolean setExtractSamplesFromOutputFiles(boolean value) {
        return configValues.put(COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, value.toString(), "boolean")
    }

    public boolean getRunFastQC() {
        return configValues.getBoolean(COConstants.FLAG_RUN_FASTQC, true)
    }

    public boolean getUseCombinedAlignAndSampe() {
        return configValues.getBoolean(COConstants.FLAG_USE_COMBINED_ALIGN_AND_SAMPE, false)
    }

    public boolean getUseExistingPairedBams() {
        return configValues.getBoolean(COConstants.FLAG_USE_EXISTING_PAIRED_BAMS, false)
    }

}
