/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.common

import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.io.ExecutionService
import groovy.transform.CompileStatic

@CompileStatic
class AlignmentAndQCConfig extends COConfig {

    public static final String CVALUE_INDEX_PREFIX = "INDEX_PREFIX"
    public static final String CVALUE_CHROMOSOME_SIZES_FILE = "CHROM_SIZES_FILE"
    public static final String CVALUE_CHROMOSOME_LENGTH_FILE = "CHROMOSOME_LENGTH_FILE_ALN"
    public static final String CVALUE_TARGET_REGIONS_FILE = "TARGET_REGIONS_FILE"
    public static final String CVALUE_TARGETSIZE = "TARGETSIZE"
    public static final String CVALUE_TARGET_SIZE = "TARGET_SIZE"
    public static final String CVALUE_RUN_COVERAGE_PLOTS_ONLY = "runCoveragePlotsOnly"
    public static final String CVALUE_OVERRIDE_BAM_FILES = "overrideBamFiles"
    public static final String CVALUE_OVERRIDE_SAMPLE_NAMES = "overrideSampleNames"
    public static final String CVALUE_MAPPABILITY_FILE = "MAPPABILITY_FILE_ALN"
    public static final String CVALUE_REPLICATION_TIME_FILE = "REPLICATION_TIME_FILE_ALN"
    public static final String CVALUE_GC_CONTENT_FILE = "GC_CONTENT_FILE_ALN"
    public static final String CVALUE_RUN_ACESEQ_QC = "runACEseqQc"

    public static final String CVALUE_CLIP_INDEX = "CLIP_INDEX"
    public static final String CVALUE_CYTOSINE_POSITIONS_INDEX = "CYTOSINE_POSITIONS_INDEX"
    public static final String CVALUE_RUN_FINGERPRINTING = "runFingerprinting"
    public static final String CVALUE_FINGERPRINTING_SITES_FILE="fingerprintingSitesFile"
    public static final String CVALUE_DEFAULT_LIBRARY_NAME = "defaultLibraryName"

    /**
     * Tool entries
     */
    public static final String TOOL_ALIGNMENT = "alignment"
    public static final String TOOL_ACCELERATED_ALIGNMENT = "accelerated:alignment"
    public static final String TOOL_COLLECT_BAM_METRICS = "collectBamMetrics"
    public static final String TOOL_SAMPESORT = "sampesort"
    public static final String TOOL_ALIGNANDPAIR = "alignAndPair"
    public static final String TOOL_ALIGNANDPAIR_SLIM = "alignAndPairSlim"
    public static final String TOOL_ACCELERATED_ALIGNANDPAIR = "accelerated:alignAndPair"
    public static final String TOOL_ACCELERATED_ALIGNANDPAIR_SLIM = "accelerated:alignAndPairSlim"
    public static final String TOOL_SAMTOOLS_INDEX = "samtoolsIndex"
    public static final String TOOL_SAMTOOLS_FLAGSTAT = "samtoolsFlagstat"
    public static final String TOOL_PURITY_ESTIMATION = "purityEstimation"
    public static final String TARGET_EXTRACTION_AND_COVERAGE_SLIM = "targetExtractCoverageSlim"

    public static final String FLAG_USE_ACCELERATED_HARDWARE = "useAcceleratedHardware"
    public static final String FLAG_USE_BIOBAMBAM_SORT = "useBioBamBamSort"
    public static final String FLAG_USE_BIOBAMBAM_MARK_DUPLICATES = "useBioBamBamMarkDuplicates"
    public static final String FLAG_USE_ONLY_EXISTING_PAIRED_BAMS = "useExistingPairedBams"
    public static final String FLAG_RUN_FASTQC = "runFastQC"
    public static final String FLAG_RUN_FASTQC_ONLY = "runFastQCOnly"
    public static final String FLAG_RUN_ALIGNMENT_ONLY = "runAlignmentOnly"
    public static final String FLAG_RUN_COVERAGE_PLOTS = "runCoveragePlots"
    public static final String FLAG_RUN_EXOME_ANALYSIS = "runExomeAnalysis"
    public static final String FLAG_RUN_COLLECT_BAMFILE_METRICS = "runCollectBamFileMetrics"
    public static final String FLAG_USE_COMBINED_ALIGN_AND_SAMPE = "useCombinedAlignAndSampe"
    public static final String FLAG_USE_SINGLE_END_PROCESSING = "useSingleEndProcessing"
    public static final String FLAG_USE_ADAPTOR_TRIMMING = "useAdaptorTrimming"
    public static final String CVALUE_MARK_DUPLICATES_VARIANT = "markDuplicatesVariant"

    AlignmentAndQCConfig(ExecutionContext context) {
        super(context)
    }

    protected String substituteDirExecution(String filename) {
        return filename.replaceAll("[\$]${ExecutionService.RODDY_CVALUE_DIRECTORY_EXECUTION}(?=\\W|\$)", dirExecution)
    }

    String getSingleBamParameter() {
        return configValues.get("bam", "")
    }

    boolean getUseOnlyExistingTargetBam() {
        return configValues.getBoolean("useOnlyExistingTargetBam", false)
    }

    boolean getUseOnlyExistingLaneBams() {
        return configValues.getBoolean(FLAG_USE_ONLY_EXISTING_PAIRED_BAMS, false)
    }

    String getIndexPrefix() {
        return configValues.getString(CVALUE_INDEX_PREFIX, "")
    }

    File getChromosomeSizesFile() {
        return new File (configValues.getString(CVALUE_CHROMOSOME_SIZES_FILE, ""))
    }

    File getChromosomeLengthFile () {
        return new File (configValues.getString(CVALUE_CHROMOSOME_LENGTH_FILE, ""))
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
        return configValues.getBoolean(FLAG_RUN_EXOME_ANALYSIS)
    }

    String getWindowSize() {
        return configValues.getString(COConstants.CVALUE_WINDOW_SIZE, "1")
    }

    boolean getRunCoveragePlotsOnly() {
        return configValues.getBoolean(CVALUE_RUN_COVERAGE_PLOTS_ONLY, false)
    }

    List<String> getOverrideMergedBamFiles() {
        return configValues.getString(CVALUE_OVERRIDE_BAM_FILES, "").split(StringConstants.SPLIT_SEMICOLON) as List<String>
    }

    List<String> getOverrideSampleNames() {
        return configValues.getString(CVALUE_OVERRIDE_SAMPLE_NAMES, "").split(StringConstants.SPLIT_SEMICOLON) as List<String>
    }

    File getMappabilityFile () {
        return new File(configValues.getString(CVALUE_MAPPABILITY_FILE))
    }

    File getReplicationTimeFile () {
        return new File(configValues.getString(CVALUE_REPLICATION_TIME_FILE))
    }

    File getGcContentFile () {
        return new File(configValues.getString(CVALUE_GC_CONTENT_FILE))
    }

    boolean getRunACEseqQC () {
        return configValues.getBoolean(CVALUE_RUN_ACESEQ_QC, false)
    }

    String getControlBamName() {
        return overrideMergedBamFiles[0]
    }

    List<String> getTumorBamNames() {
        return overrideMergedBamFiles[1..-1]
    }

    void setSampleExtractionFromOutputFiles(boolean value) {
        setConfig(COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, value.toString(), "boolean")
    }

    File getCytosinePositionIndex() {
        return new File(configValues.getString(CVALUE_CYTOSINE_POSITIONS_INDEX))
    }

    boolean getUseAdapterTrimming() {
        return configValues.getBoolean(FLAG_USE_ADAPTOR_TRIMMING, false)
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

    String getDirExecution() {
        return context.executionDirectory.absolutePath
    }

    boolean getRunFastQCOnly() {
        return configValues.getBoolean(FLAG_RUN_FASTQC_ONLY, false)
    }

    boolean getRunAlignmentOnly() {
        return configValues.getBoolean(FLAG_RUN_ALIGNMENT_ONLY, false)
    }

    boolean getRunCoveragePlots() {
        return configValues.getBoolean(FLAG_RUN_COVERAGE_PLOTS, true)
    }

    boolean getRunCollectBamFileMetrics() {
        return configValues.getBoolean(FLAG_RUN_COLLECT_BAMFILE_METRICS, false)
    }

    boolean setExtractSamplesFromOutputFiles(boolean value) {
        return configValues.put(COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, value.toString(), "boolean")
    }

    boolean getRunFastQC() {
        return configValues.getBoolean(FLAG_RUN_FASTQC, true)
    }

    boolean getUseCombinedAlignAndSampe() {
        return configValues.getBoolean(FLAG_USE_COMBINED_ALIGN_AND_SAMPE, false)
    }

    String getDefaultLibraryName() {
        return configValues.getString(CVALUE_DEFAULT_LIBRARY_NAME, "null")
    }

    String getUseSingleEndProcessing() {
        return configValues.getBoolean(FLAG_USE_SINGLE_END_PROCESSING , false)
    }

}
