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
    public static final String CVALUE_CHROMOSOME_LENGTH_FILE = "CHROMOSOME_LENGTH_FILE"
    public static final String CVALUE_TARGET_REGIONS_FILE = "TARGET_REGIONS_FILE"
    public static final String CVALUE_TARGETSIZE = "TARGETSIZE"
    public static final String CVALUE_TARGET_SIZE = "TARGET_SIZE"
    public static final String CVALUE_RUN_COVERAGE_PLOTS_ONLY = "runCoveragePlotsOnly";
    public static final String CVALUE_OVERRIDE_BAM_FILES = "overrideBamFiles";
    public static final String CVALUE_OVERRIDE_SAMPLE_NAMES = "overrideSampleNames";
    public static final String CVALUE_MAPPABILITY_FILE = "MAPPABILITY_FILE";
    public static final String CVALUE_REPLICATION_TIME_FILE = "REPLICATION_TIME_FILE";
    public static final String CVALUE_GC_CONTENT_FILE = "GC_CONTENT_FILE";
    public static final String CVALUE_RUN_ACESEQ_QC = "runACEseqQc";

    public static final String CVALUE_CLIP_INDEX = "CLIP_INDEX"
    public static final String CVALUE_CYTOSINE_POSITIONS_INDEX = "CYTOSINE_POSITIONS_INDEX"

    public AlignmentAndQCConfig(ExecutionContext context) {
        super(context)
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

    public File getChromosomeLengthFile () {
        return new File (configValues.getString(CVALUE_CHROMOSOME_LENGTH_FILE, ""))
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

    public boolean getRunFastQCOnly() {
        return configValues.getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false)
    }

    public boolean getRunFastqQC() {
        return configValues.getBoolean(COConstants.FLAG_RUN_FASTQC, true)
    }

    public boolean getRunAlignmentOnly() {
        return configValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false)
    }

    public boolean getRunCoveragePlots() {
        return configValues.getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true)
    }

    public boolean getRunCollectBamFileMetrics() {
        return configValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false)
    }

    public boolean getUseExistingPairedBams() {
        return configValues.getBoolean(COConstants.FLAG_USE_EXISTING_PAIRED_BAMS, false)
    }

    public boolean getUseCombinedAlignAndSampe() {
        return configValues.getBoolean(COConstants.FLAG_USE_COMBINED_ALIGN_AND_SAMPE, false)
    }

    public boolean getRunSlimWorkflow() {
        return configValues.getBoolean(COConstants.FLAG_RUN_SLIM_WORKFLOW, false)
    }

    public String getWindowSize() {
        return configValues.getString(COConstants.CVALUE_WINDOW_SIZE, "1")
    }

    public boolean getRunCoveragePlotsOnly() {
        return configValues.getBoolean(CVALUE_RUN_COVERAGE_PLOTS_ONLY, false)
    }

    public List<String> getOverrideMergedBamFiles() {
        return configValues.getString(CVALUE_OVERRIDE_BAM_FILES, "").split(StringConstants.SPLIT_SEMICOLON) as List<String>
    }

    public List<String> getOverrideSampleNames() {
        return configValues.getString(CVALUE_OVERRIDE_SAMPLE_NAMES, "").split(StringConstants.SPLIT_SEMICOLON) as List<String>
    }

    public File getMappabilityFile () {
        return new File(configValues.getString(CVALUE_MAPPABILITY_FILE))
    }

    public File getReplicationTimeFile () {
        return new File(configValues.getString(CVALUE_REPLICATION_TIME_FILE))
    }

    public File getGcContentFile () {
        return new File(configValues.getString(CVALUE_GC_CONTENT_FILE))
    }

    public boolean getRunACEseqQC () {
        return configValues.getBoolean(CVALUE_RUN_ACESEQ_QC, false)
    }

    public String getControlBamName() {
        return overrideMergedBamFiles[0]
    }

    public List<String> getTumorBamNames() {
        return overrideMergedBamFiles[1..-1]
    }

    public void setSampleExtractionFromOutputFiles(boolean value) {
        setConfig(COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, value.toString(), "boolean")
    }

    public File getCytosinePositionIndex() {
        return new File(configValues.getString(CVALUE_CYTOSINE_POSITIONS_INDEX))
    }

    public File getClipIndex() {
        return new File(configValues.getString(CVALUE_CLIP_INDEX))
    }

}
