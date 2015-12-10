package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.files.*;
import de.dkfz.b080.co.common.*
import de.dkfz.roddy.core.*
import de.dkfz.roddy.knowledge.files.Tuple3

import java.util.*

/**
 * @author michael
 */
public class QCPipeline extends Workflow {

    public QCPipeline() {}

    @Override
    public boolean execute(ExecutionContext context) {
        QCConfig cfg = new QCConfig(context)
        cfg.sampleExtractionFromOutputFiles = false
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getProject().getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForContext(context);
        if (samples.size() == 0)
            return false;

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>();

        for (Sample sample : samples) {
            BamFileGroup sortedBamFiles = createSortedBams(cfg, runtimeService, sample);

            if (sortedBamFiles.getFilesInGroup().size() == 0) continue;

            if (cfg.runFastQCOnly || cfg.runAlignmentOnly) continue;

            BamFile mergedBam;
            if (cfg.runSlimWorkflow) {
                mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample);
                if (cfg.runCollectBamFileMetrics) mergedBam.collectMetrics();
            } else {
                mergedBam = mergeAndRemoveDuplicatesFat(cfg, sample, sortedBamFiles);
            }

            if (cfg.runExomeAnalysis) {
                if(!cfg.runSlimWorkflow) mergedBam.rawBamCoverage();
                BamFile targetOnlyBamFile = mergedBam.extractTargetsCalculateCoverage();
            }

            Sample.SampleType sampleType = sample.getType();
            if (!coverageTextFilesBySample.containsKey(sampleType))
                coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup());
            coverageTextFilesBySample.get(sampleType).addFile(mergedBam.calcReadBinsCoverage());

            mergedBamFiles.addFile(mergedBam);
        }

        if (mergedBamFiles.getFilesInGroup().size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."));
            return false;
        }

        if (cfg.runFastQCOnly)
            return true;

        if (cfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot();
        }

        return true;
    }

    private boolean runSlim(ExecutionContext context) {
        return false;
    }

    private boolean runFat(ExecutionContext context) {
        return false;
    }

    /**
     * This entry is used for caching purposes.
     */
    protected Map<DataSet, Map<String, List<LaneFileGroup>>> foundRawSequenceFileGroups = new LinkedHashMap<>();

    /**
     * Provides a cached method for loading lane files from a sample.
     *
     * @param sample
     * @return
     */
    protected synchronized List<LaneFileGroup> loadLaneFilesForSample(ExecutionContext context, Sample sample) {
        DataSet dataSet = context.getDataSet();
        if (!foundRawSequenceFileGroups.containsKey(dataSet)) {
            foundRawSequenceFileGroups.put(dataSet, new LinkedHashMap<String, List<LaneFileGroup>>());
        }
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService();
        String sampleID = sample.getName();
        Map<String, List<LaneFileGroup>> mapForDataSet = foundRawSequenceFileGroups.get(dataSet);
        if (!mapForDataSet.containsKey(sampleID)) {
            List<LaneFileGroup> laneFileGroups = runtimeService.getLanesForSample(context, sample);
            mapForDataSet.put(sampleID, laneFileGroups);
        }

        List<LaneFileGroup> laneFileGroups = mapForDataSet.get(sampleID);
        List<LaneFileGroup> copyOfLaneFileGroups = new LinkedList<LaneFileGroup>();
        for (LaneFileGroup lfg : laneFileGroups) {
            List<LaneFile> copyOfFiles = new LinkedList<>();
            for (LaneFile lf : lfg.getFilesInGroup()) {
                LaneFile copyOfFile = new LaneFile(lf, context);//lf.getPath(), context, lf.getCreatingJobsResult(), lf.getParentFiles(), lf.getFileStage());
                copyOfFiles.add(copyOfFile);
            }
            copyOfLaneFileGroups.add(new LaneFileGroup(context, lfg.getId(), lfg.getRun(), sample, copyOfFiles));
        }
        return copyOfLaneFileGroups;
    }

    private void aceSeqQc(CoverageTextFile windowedCoverageTextFile) {
        TextFile annotationResult = AceSeqQC.annotateCovWindows(windowedCoverageTextFile);
        TextFile mergedAndFilteredCoverageWindowFiles = AceSeqQC.mergeAndFilterCovWindows(annotationResult);
        Tuple3<TextFile, TextFile, TextFile> correctedWindowFile = AceSeqQC.correctGC(mergedAndFilteredCoverageWindowFiles);
    }

    private BamFileGroup createSortedBams(QCConfig cfg, COProjectsRuntimeService runtimeService, Sample sample) {
        BamFileGroup sortedBamFiles = new BamFileGroup();

        if (cfg.useExistingPairedBams) {
            //Start from the paired bams instead of the lane files.
            sortedBamFiles = runtimeService.getPairedBamFilesForDataSet(cfg.context, sample);

        } else {
            //Create bam files out of the lane files
            List<LaneFileGroup> rawSequenceGroups = loadLaneFilesForSample(cfg.context, sample);
            if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
                return sortedBamFiles;
            for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                if (cfg.runFastqQC && !cfg.runAlignmentOnly)
                    rawSequenceGroup.calcFastqcForAll();
                if (cfg.runFastQCOnly)
                    continue;

                BamFile bamFile = null;

                if (cfg.useCombinedAlignAndSampe) { //I.e. bwa mem
                    if (cfg.runSlimWorkflow) {
                        bamFile = rawSequenceGroup.alignAndPairSlim();
                    } else {
                        bamFile = rawSequenceGroup.alignAndPair();
                    }
                } else { //I.e. bwa align
                    rawSequenceGroup.alignAll();
                    if (cfg.runSlimWorkflow) {
                        bamFile = rawSequenceGroup.getAllAlignedFiles().pairAndSortSlim();
                    } else {
                        bamFile = rawSequenceGroup.getAllAlignedFiles().pairAndSort();
                    }
                }

                // @Michael: The comments suggests that this should only be called in the else branch above!?
                // @Michael: Why should BAM files created with sai files be temporary?
                bamFile.setAsTemporaryFile();  // Bam files created with sai files are only temporary.
                sortedBamFiles.addFile(bamFile);

                if (cfg.windowSize == 1) {
                    aceSeqQc(bamFile.readBinsCoverageTextFile)
                } else {
                    // throw "This won't work with s.th. else than 1kb!
                    // TODO commonCOWorkflowSettings: 10kb, exome 10kb
                }
            }

        }
        return sortedBamFiles;
    }

    private BamFile mergeAndRemoveDuplicatesFat(QCConfig cfg, Sample sample, BamFileGroup sortedBamFiles) {

        //To avoid problems with qcsummary the step is done manually.
        sortedBamFiles.runDefaultOperations();
        for (BamFile sortedBam : sortedBamFiles.getFilesInGroup())
            sortedBam.createQCSummaryFile();

        BamFile mergedBam = sortedBamFiles.mergeAndRemoveDuplicates();

        if (cfg.runCollectBamFileMetrics) mergedBam.collectMetrics();
        mergedBam.runDefaultOperations();
        mergedBam.calcCoverage();

        Sample.SampleType sampleType = sample.getType();


        mergedBam.createQCSummaryFile();
        return mergedBam;
    }

    @Override
    public boolean checkExecutability(ExecutionContext context) {
        QCConfig cfg = new QCConfig(context)
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getProject().getRuntimeService();
        List<Sample> samples = runtimeService.getSamplesForContext(context);
        if (samples.size() == 0)
            return false;

        if (!cfg.useExistingPairedBams) {
            //Check if at least one file is available. Maybe for two if paired is used...?
            int cnt = 0;
            for (Sample sample : samples) {

                List<LaneFileGroup> laneFileGroups = loadLaneFilesForSample(context, sample);
                for (LaneFileGroup lfg : laneFileGroups) {
                    cnt += lfg.getFilesInGroup().size();
                }
            }
            return cnt > 0;
        } else {
            return true;
        }
    }

    @Override
    public boolean createTestdata(ExecutionContext context) {
        boolean allOk = true;
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getProject().getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForContext(context);
        for (Sample sample : samples) {
            List<LaneFile> files = new LinkedList<LaneFile>();
            LaneFileGroup allLaneFiles = new LaneFileGroup(context, "allLaneFiles", "noSpecificRun", sample, files);

            List<LaneFileGroup> rawSequenceGroups = runtimeService.getLanesForSample(context, sample);
            for (LaneFileGroup lfg : rawSequenceGroups) {
                for (LaneFile lf : lfg.getFilesInGroup()) {
                    allLaneFiles.addFile(lf);
                }
            }
            allLaneFiles.createTestDataForLaneFiles();
        }
        return allOk;
    }
}
