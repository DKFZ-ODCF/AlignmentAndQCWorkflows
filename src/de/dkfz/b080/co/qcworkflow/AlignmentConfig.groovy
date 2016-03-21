package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.common.COConfig
import de.dkfz.b080.co.files.AlignmentConstants
import de.dkfz.roddy.core.ExecutionContext

// TODO Discuss which Constants can be moved into the alignment plugin
// TODO Discuss whether configuration value constants should be in AlignmentConfig/COConfig rather than in dedicated class.
import static de.dkfz.b080.co.files.AlignmentConstants.*
import static de.dkfz.b080.co.files.COConstants.*
import de.dkfz.roddy.StringConstants


@groovy.transform.CompileStatic
class AlignmentConfig extends COConfig {

    public AlignmentConfig (ExecutionContext context) {
        super(context)
    }

    public boolean getRunFastQCOnly() {
        return configValues.getBoolean(FLAG_RUN_FASTQC_ONLY, false)
    }

    public boolean getRunFastqQC() {
        return configValues.getBoolean(FLAG_RUN_FASTQC, true)
    }

    public boolean getRunAlignmentOnly() {
        return configValues.getBoolean(FLAG_RUN_ALIGNMENT_ONLY, false)
    }

    public boolean getRunCoveragePlots() {
        return configValues.getBoolean(FLAG_RUN_COVERAGE_PLOTS, true)
    }

    public boolean getRunExomeAnalysis() {
        return configValues.getBoolean(FLAG_RUN_EXOME_ANALYSIS)
    }

    public boolean getRunCollectBamFileMetrics() {
        return configValues.getBoolean(FLAG_RUN_COLLECT_BAMFILE_METRICS, false)
    }

    public boolean getUseExistingPairedBams() {
        return configValues.getBoolean(FLAG_USE_EXISTING_PAIRED_BAMS, false)
    }

    public boolean getUseCombinedAlignAndSampe() {
        return configValues.getBoolean(FLAG_USE_COMBINED_ALIGN_AND_SAMPE, false)
    }

    public boolean getRunSlimWorkflow() {
        return configValues.getBoolean(FLAG_RUN_SLIM_WORKFLOW, false)
    }

    public String getWindowSize() {
        return configValues.getString(CVALUE_WINDOW_SIZE, "1")
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

    public String getMappabilityFile () {
        return configValues.getString(CVALUE_MAPPABILITY_FILE)
    }

    public String getReplicationTimeFile () {
        return configValues.getString(CVALUE_REPLICATION_TIME_FILE)
    }

    public String getGcContentFile () {
        return configValues.getString(CVALUE_GC_CONTENT_FILE)
    }

    public boolean getRunACEseqQC () {
        return configValues.getBoolean(CVALUE_RUN_ACESEQ_QC, true)
    }

    public String getControlBamName() {
        return overrideMergedBamFiles[0]
    }

    public List<String> getTumorBamNames() {
        return overrideMergedBamFiles[1..-1]
    }

    public void setSampleExtractionFromOutputFiles(boolean value) {
        setConfig(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, value.toString(), "boolean" )
    }

}
