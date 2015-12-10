package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.files.COConstants
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.core.ExecutionContext

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES
import static de.dkfz.b080.co.files.COConstants.FLAG_IS_DEBUG_CONFIGURATION;

class QCConfig {

    public final ExecutionContext context

    QCConfig(ExecutionContext context) {
        this.context = context
    }

    public setConfig(String flagName, String value, String typeName) {
        extractConfigValue(context).put(flagName, value, typeName);
    }

    // This is used so often, it should be part of ExecutionContext.
    public extractConfigValue(ExecutionContext context) {
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
        extractConfigValue(context).getBoolean("runCoveragePlotsOnly", false)
    }

    public boolean getRunSlimWorkflow() {
        extractConfigValue(context).getBoolean(COConstants.FLAG_RUN_SLIM_WORKFLOW, false)
    }

    public String[] getWindowSize() {
        extractConfigValue(context).getString(COConstants.CVALUE_WINDOW_SIZE, 1)
    }

    public String getOverrideFastqFiles() {
        extractConfigValue(context).getString("overrideFastqFiles", "")
    }

    public String[] getOverrideMergedBamFiles() {
        extractConfigValue(context).getString("overrideBamFiles", "").split(StringConstants.SPLIT_SEMICOLON)
    }

    public String[] getOverrideSampleNames() {
        extractConfigValue(context).getString("overrideSampleNames", "").split(StringConstants.SPLIT_SEMICOLON)
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



}
