package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.files.AlignmentConstants
import de.dkfz.b080.co.files.COConstants
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.core.ExecutionContext

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES

class QCConfig {

    public final ExecutionContext context

    QCConfig(ExecutionContext context) {
        this.context = context
    }

    public setConfig(String flagName, String value, String typeName) {
        extractConfigValue(context).put(flagName, value, typeName);
    }

    // This is used so often, it should be part of ExecutionContext.
    public Object extractConfigValue(ExecutionContext context) {
        context.getConfiguration().getConfigurationValues();
    }

    public boolean getRunFastQCOnly() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false)
    }

    public boolean getRunFastqQC() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_FASTQC, true)
    }

    public boolean getRunAlignmentOnly() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false)
    }

    public boolean getRunCoveragePlots() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true)
    }

    public boolean getRunExomeAnalysis() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_EXOME_ANALYSIS)
    }

    public boolean getRunCollectBamFileMetrics() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false)
    }

    public boolean getUseExistingPairedBams() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_USE_EXISTING_PAIRED_BAMS, false)
    }

    public boolean getUseCombinedAlignAndSampe() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_USE_COMBINED_ALIGN_AND_SAMPE, false)
    }

    public boolean getRunCoveragePlotsOnly() {
        extractConfigValue(context).getBoolean(AlignmentConstants.CVALUE_RUN_COVERAGE_PLOTS_ONLY, false)
    }

    public boolean getRunSlimWorkflow() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_SLIM_WORKFLOW, false)
    }

    public String getWindowSize() {
        extractConfigValue(context).getString(COConstants.CVALUE_WINDOW_SIZE, "1")
    }

    public String[] getOverrideFastqFiles() {
        extractConfigValue(context).getString(AlignmentConstants.CVALUE_OVERRIDE_FASTQ_FILES, "")
    }

    public String[] getOverrideMergedBamFiles() {
        extractConfigValue(context).getString(AlignmentConstants.CVALUE_OVERRIDE_BAM_FILES, "").split(StringConstants.SPLIT_SEMICOLON)
    }

    public String[] getOverrideSampleNames() {
        extractConfigValue(context).getString(AlignmentConstants.CVALUE_OVERRIDE_SAMPLE_NAMES, "").split(StringConstants.SPLIT_SEMICOLON)
    }

    public String getControlBamName() {
        overrideMergedBamFiles[0]
    }

    public String[] getTumorBamNames() {
        overrideMergedBamFiles[1..-1]
    }

    public void setSampleExtractionFromOutputFiles(boolean value) {
        setConfig(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, value.toString(), "boolean" )
    }

    public String getMappabilityFile () {
        extractConfigValue(context).getString(AlignmentConstants.CVALUE_MAPPABILITY_FILE)
    }

    public String getReplicationTimeFile () {
        extractConfigValue(context).getString(AlignmentConstants.CVALUE_REPLICATION_TIME_FILE)
    }

    public String getGcContentFile () {
        extractConfigValue(context).getString(AlignmentConstants.CVALUE_GC_CONTENT_FILE)
    }

    public boolean getRunACEseqQC () {
        extractConfigValue(context).getBoolean(AlignmentConstants.CVALUE_RUN_ACESEQ_QC, true)
    }

}
