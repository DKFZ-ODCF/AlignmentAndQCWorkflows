/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

import de.dkfz.b080.co.common.AlignmentAndQCConfig;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.jobs.Job;
import de.dkfz.roddy.execution.jobs.ScriptCallingMethod;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.ITestdataSource;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * A bam file is the binary version of a sam file and contains sequence data.
 *
 * @author heinold
 * @url http://www.broadinstitute.org/igv/bam
 */
public class BamFile extends BasicBamFile implements ITestdataSource {

    private BamIndexFile indexFile;
    private FlagstatsFile flagstatsFile;
    private TextFile extendedFlagstatsFile;
    private ChromosomeDiffFileGroup chromosomeDiffFileGroup;
    private ChromosomeDiffValueFile chromosomeDiffStatisticsFile;
    private ChromosomeDiffTextFile chromosomeDiffMatrixFile;
    private ChromosomeDiffPlotFile chromosomeDiffPlotFile;
    private InsertSizesFileGroup insertSizesFileGroup;
    private InsertSizesValueFile insertSizesStatisticsFile;
    private InsertSizesTextFile insertSizesMatrixFile;
    private InsertSizesPlotFile insertSizesPlotFile;
    private CoverageTextFile genomeCoverageTextFile;
    private CoverageTextFile targetCoverageTextFile;
    private CoverageTextFile rawBamCoverageTextFile;
    private CoverageTextFile readBinsCoverageTextFile;
    private CoverageTextFile targetsWithCoverageTextFile;
    private TextFile groupedGenomeCoverageTextFile;
    private BamMetricsFile bamMetricsFile;
    private BamFile targetExtractedBamFile;
    private boolean _isTargetExtractedFile;
    private OnTargetCoverageTextFile onTargetCoverageTextFile;
    private QCSummaryFile qcSummaryFile;
    private TextFile qcJsonFile;
    private MethylationMetaCheckpointFile methylationMetaCheckpointFile;
    private MethylationMetaMetricsCheckpointFile methylationMetaMetricsCheckpointFile;
    private TextFile fingerprintsFile;


    public BamFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    /**
     * "Copy" constructor for Bam conversion.
     *
     * @param minorBam
     */
    public BamFile(BasicBamFile minorBam) {
        super(minorBam);
    }

    public String getType() {
        return _isTargetExtractedFile ? "targetExtractedBam" : "commonBam";
    }

    public Sample getSample() {
        return ((COFileStageSettings) getFileStage()).getSample();
    }

    public boolean hasIndex() {
        return indexFile != null;
    }

    public boolean hasFlagstats() {
        return flagstatsFile != null;
    }

    public boolean hasInsertSizes() {
        return insertSizesFileGroup != null;
    }

    public boolean hasChromosomeDiff() {
        return chromosomeDiffFileGroup != null;
    }

    public boolean hasCoverage() {
        return genomeCoverageTextFile != null;
    }

    public boolean hasMetricsFile() {
        return bamMetricsFile != null;
    }

    public void setMetricsFile(BamMetricsFile metricsFile) {
        this.bamMetricsFile = metricsFile;
    }

    public boolean hasOnTargetCoverageTextFile() {
        return onTargetCoverageTextFile != null;
    }

    public boolean isRawBamFile() {
        return !_isTargetExtractedFile;
    }

    public void setIsTargetExtracted() {
        _isTargetExtractedFile = true;
    }

    public boolean isTargetExtractedBamFile() {
        return _isTargetExtractedFile;
    }

    public boolean hasTargetExtractedBamFile() {
        return targetExtractedBamFile != null;
    }

    public boolean hasRawBamCoverageTextFile() {
        return rawBamCoverageTextFile != null;
    }

    public boolean hasTargetCoverageTextFile() {
        return targetCoverageTextFile != null;
    }

    public void setIndexFile(BamIndexFile indexFile) {
        this.indexFile = indexFile;
    }

    public void setFlagstatsFile(FlagstatsFile flagstatsFile) {
        this.flagstatsFile = flagstatsFile;
    }

    public BamIndexFile getIndexFile() {
        return indexFile;
    }

    public ChromosomeDiffFileGroup getChromosomeDiffFileGroup() {
        return chromosomeDiffFileGroup;
    }

    public CoverageTextFile getGenomeCoverageTextFile() { return genomeCoverageTextFile; }

    public TextFile getGroupedGenomeCoverageTextFile() { return groupedGenomeCoverageTextFile; }

    public FlagstatsFile getFlagstatsFile() {
        return flagstatsFile;
    }

    public void setExtendedFlagstatsFile(TextFile extendedFlagstatsFile) {
        this.extendedFlagstatsFile = extendedFlagstatsFile;
    }

    public TextFile getExtendedFlagstatsFile() {
        return extendedFlagstatsFile;
    }

    public BamMetricsFile getMetricsFile() {
        return bamMetricsFile;
    }

    public InsertSizesFileGroup getInsertSizesFileGroup() {
        return insertSizesFileGroup;
    }

    public CoverageTextFile getTargetCoverageTextFile() {
        return targetCoverageTextFile;
    }

    public CoverageTextFile getRawBamCoverageTextFile() {
        return rawBamCoverageTextFile;
    }

    public MethylationMetaCheckpointFile getMethylationMetaCheckpointFile() {
        return methylationMetaCheckpointFile;
    }

    public MethylationMetaMetricsCheckpointFile getMethylationMetaMetricsCheckpointFile() { return methylationMetaMetricsCheckpointFile; }

    public QCSummaryFile getQcSummaryFile() {
        return qcSummaryFile;
    }

    public TextFile getQcJsonFile() {
        return qcJsonFile;
    }

    public void setChromosomeDiffStatisticsFile(ChromosomeDiffValueFile chromosomeDiffStatisticsFile) {
        this.chromosomeDiffStatisticsFile = chromosomeDiffStatisticsFile;
    }

    public void setChromosomeDiffMatrixFile(ChromosomeDiffTextFile chromosomeDiffMatrixFile) {
        this.chromosomeDiffMatrixFile = chromosomeDiffMatrixFile;
    }

    public void setChromosomeDiffPlotFile(ChromosomeDiffPlotFile chromosomeDiffPlotFile) {
        this.chromosomeDiffPlotFile = chromosomeDiffPlotFile;
    }

    public void setInsertSizesStatisticsFile(InsertSizesValueFile insertSizesStatisticsFile) {
        this.insertSizesStatisticsFile = insertSizesStatisticsFile;
    }

    public void setInsertSizesMatrixFile(InsertSizesTextFile insertSizesMatrixFile) {
        this.insertSizesMatrixFile = insertSizesMatrixFile;
    }

    public void setInsertSizesPlotFile(InsertSizesPlotFile insertSizesPlotFile) {
        this.insertSizesPlotFile = insertSizesPlotFile;
    }

    public void setGenomeCoverageTextFile(CoverageTextFile genomeCoverageTextFile) {
        this.genomeCoverageTextFile = genomeCoverageTextFile;
    }

    public void setGroupedGenomeCoverageTextFile(TextFile groupedGenomeCoverageTextFile) { this.groupedGenomeCoverageTextFile = groupedGenomeCoverageTextFile; }

    public void setReadBinsCoverageTextFile(CoverageTextFile readBinsCoverageTextFile) {
        this.readBinsCoverageTextFile = readBinsCoverageTextFile;
    }

    public void setBamMetricsFile(BamMetricsFile bamMetricsFile) {
        this.bamMetricsFile = bamMetricsFile;
    }

    public void setQcSummaryFile(QCSummaryFile qcSummaryFile) {
        this.qcSummaryFile = qcSummaryFile;
    }

    public void setQcJsonFile(TextFile qcJsonFile) {
        this.qcJsonFile = qcJsonFile;
    }

    public ChromosomeDiffValueFile getChromosomeDiffStatisticsFile() {
        return chromosomeDiffStatisticsFile;
    }

    public ChromosomeDiffTextFile getChromosomeDiffMatrixFile() {
        return chromosomeDiffMatrixFile;
    }

    public ChromosomeDiffPlotFile getChromosomeDiffPlotFile() {
        return chromosomeDiffPlotFile;
    }

    public InsertSizesValueFile getInsertSizesStatisticsFile() {
        return insertSizesStatisticsFile;
    }

    public InsertSizesTextFile getInsertSizesMatrixFile() {
        return insertSizesMatrixFile;
    }

    public InsertSizesPlotFile getInsertSizesPlotFile() {
        return insertSizesPlotFile;
    }

    public BamMetricsFile getBamMetricsFile() {
        return bamMetricsFile;
    }

    public BamFile getTargetExtractedBamFile() {
        return targetExtractedBamFile;
    }

    public TextFile getFingerprintsFile() { return fingerprintsFile; }

    public void setTargetCoverageTextFile(CoverageTextFile targetCoverageTextFile) {
        this.targetCoverageTextFile = targetCoverageTextFile;
    }

    public CoverageTextFile getTargetsWithCoverageTextFile() {
        return targetsWithCoverageTextFile;
    }

    public void setFingerprintsFile(TextFile fingerprintsFile) {
        this.fingerprintsFile = fingerprintsFile;
    }

    public void setTargetsWithCoverageTextFile(CoverageTextFile targetsWithCoverageTextFile) {
        this.targetsWithCoverageTextFile = targetsWithCoverageTextFile;
    }

    public void setOnTargetCoverageTextFile(OnTargetCoverageTextFile tf) {
        this.onTargetCoverageTextFile = tf;
    }

    public OnTargetCoverageTextFile getOnTargetCoverageTextFile() {
        return onTargetCoverageTextFile;
    }

    @ScriptCallingMethod
    public Tuple2<MethylationMetaMetricsCheckpointFile,MethylationMetaCheckpointFile> methylationCallingMeta() {
        if (methylationMetaCheckpointFile == null || methylationMetaMetricsCheckpointFile == null) {
            Tuple2<MethylationMetaMetricsCheckpointFile, MethylationMetaCheckpointFile> res = GenericMethod.callGenericTool("methylationCallingMeta", this);
            methylationMetaCheckpointFile = res.value1;
            methylationMetaMetricsCheckpointFile = res.value0;
        }
        return new Tuple2(methylationMetaCheckpointFile, methylationMetaMetricsCheckpointFile);
    }

    @ScriptCallingMethod
    public Tuple2<MethylationMetaMetricsCheckpointFile,MethylationMetaCheckpointFile> libraryMethylationCallingMeta() {
        return methylationCallingMeta();
    }

    @ScriptCallingMethod
    public Tuple2<MethylationMetaMetricsCheckpointFile,MethylationMetaCheckpointFile>  mergedMethylationCallingMeta() {
        return methylationCallingMeta();
    }

    // This code has not been tested for ages.
    @Deprecated
    @ScriptCallingMethod
    public BamMetricsCollectionFileGroup collectMetrics() {
        final ExecutionContext run = getExecutionContext();
        final Configuration cfg = run.getConfiguration();
        BamMetricsCollectionFileGroup resultGroup = new BamMetricsCollectionFileGroup(new LinkedList());

        String metricsOutputDirectory = cfg.getConfigurationValues().getString("metricsOutputDirectory");

        List<BaseFile> pFiles = Arrays.asList((BaseFile) this);

        final String TOOL = AlignmentAndQCConfig.TOOL_COLLECT_BAM_METRICS;
        Map<String, Object> parameters = run.getDefaultJobParameters(TOOL);
        parameters.put("FILENAME", this.getPath().getAbsolutePath());
        parameters.put("DIR_METRICS", run.getRuntimeService().getDirectory(metricsOutputDirectory, run).getAbsolutePath());
        parameters.put("sample", ((COFileStageSettings) this.getFileStage()).getSample().getName());
        parameters.put("SAMPLE", ((COFileStageSettings) this.getFileStage()).getSample().getName());

        BamMetricsAlignmentSummaryFile alnSummaryFile = new BamMetricsAlignmentSummaryFile(this);
        resultGroup.addFile(alnSummaryFile);

        Job job = new Job(run, run.createJobName(this, TOOL), TOOL, parameters, pFiles, Arrays.asList((BaseFile) alnSummaryFile));
        alnSummaryFile.setCreatingJobsResult(job.run());
        return resultGroup;
    }

    @ScriptCallingMethod
    public BamFile extractTargetsCalculateCoverage() {
        if (targetExtractedBamFile == null)
            targetExtractedBamFile = GenericMethod.callGenericTool(COConstants.TARGET_EXTRACTION_AND_COVERAGE_SLIM, this, this.getGenomeCoverageTextFile(), "SAMPLE=" + getSample().getName());
        return targetExtractedBamFile;
    }

    @ScriptCallingMethod
    public CoverageTextFile getReadBinsCoverageTextFile() {
        return calcReadBinsCoverage();
    }

    /**
     * Counts reads in Window-Size kb
     */
    @ScriptCallingMethod
    public CoverageTextFile calcReadBinsCoverage() {
        if (readBinsCoverageTextFile == null)
            readBinsCoverageTextFile = GenericMethod.callGenericTool("readBinsCoverage", this);
        return readBinsCoverageTextFile;
    }

    public List<COBaseFile> getCreatedFiles() {
        List<COBaseFile> files = new LinkedList<>();
        if (hasIndex()) {
            files.add(getIndexFile());
        }
        if (hasFlagstats()) {
            files.add(flagstatsFile);
        }
        if (hasChromosomeDiff()) {
            files.add(getChromosomeDiffFileGroup().getTextFile());
            files.add(getChromosomeDiffFileGroup().getPlotFile());
            files.add(getChromosomeDiffFileGroup().getValueFile());
        }
        if (hasCoverage()) {
            files.add(getGenomeCoverageTextFile());
        }
        if (hasInsertSizes()) {
            files.add(getInsertSizesFileGroup().getTextFile());
            files.add(getInsertSizesFileGroup().getPlotFile());
            files.add(getInsertSizesFileGroup().getValueFile());
        }
        if (hasMetricsFile()) {
            files.add(getMetricsFile());
        }
        //TODO Find out for which files / operations an indexFile is used
        if (hasTargetCoverageTextFile()) {
            files.add(getTargetCoverageTextFile());
        }
        if (hasRawBamCoverageTextFile()) {
            files.add(getRawBamCoverageTextFile());
        }
        if (hasTargetExtractedBamFile()) {
            files.add(getTargetExtractedBamFile());
        }
        if (hasOnTargetCoverageTextFile()) {
            files.add(getOnTargetCoverageTextFile());
        }

        return files;
    }

    @Override
    public String toString() {
        return String.format("BamFile %s hasIndex=%b hasFlagstats=%b hasChromosomeDiff=%b hasCoverage=%b",
                super.toString(), hasIndex(), hasFlagstats(), hasChromosomeDiff(), hasCoverage());
    }

    @Override
    public boolean createTestData() {
        throw new UnsupportedOperationException("Creating test data for BamFile not supported yet.");
    }
}
