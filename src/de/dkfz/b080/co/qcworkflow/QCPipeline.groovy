package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.files.*;
import de.dkfz.b080.co.common.*
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.tools.LoggerWrapper;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues;
import de.dkfz.roddy.core.*

import java.util.*;

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES;

/**
 * @author michael
 */
@groovy.transform.CompileStatic
public class QCPipeline extends Workflow {

    private static LoggerWrapper logger = LoggerWrapper.getLogger(QCPipeline.class.getName());

    public QCPipeline() {}

    @Override
    public boolean execute(ExecutionContext context) {
        Configuration cfg = context.getConfiguration();
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context);

        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues();
        cfgValues.put(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, "false", "boolean"); //Disable sample extraction from output for alignment workflows.

        // Run flags
        final boolean runFastQCOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false);
        final boolean runAlignmentOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false);
        final boolean runCoveragePlots = cfgValues.getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true);
        final boolean runSlimWorkflow = cfgValues.getBoolean(COConstants.FLAG_RUN_SLIM_WORKFLOW, false);
        final boolean runExomeAnalysis = cfgValues.getBoolean(COConstants.FLAG_RUN_EXOME_ANALYSIS);
        final boolean runCollectBamFileMetrics = cfgValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false);

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForContext(context);
        if (samples.size() == 0)
            return false;

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>();

        for (Sample sample : samples) {
            BamFileGroup sortedBamFiles = createSortedBams(context, runtimeService, sample);

            if (sortedBamFiles.getFilesInGroup().size() == 0) continue;

            if (runFastQCOnly || runAlignmentOnly) continue;

            BamFile mergedBam;
            if (runSlimWorkflow) {
                mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample);
                if (runCollectBamFileMetrics) mergedBam.collectMetrics();
            } else {
                mergedBam = mergeAndRemoveDuplicatesFat(context, sample, sortedBamFiles);
            }

            if (runExomeAnalysis) {
                if (!runSlimWorkflow) mergedBam.rawBamCoverage();
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

        if (runFastQCOnly)
            return true;

        if (runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot();
        }

        return true;
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

    private BamFileGroup createSortedBams(ExecutionContext context, COProjectsRuntimeService runtimeService, Sample sample) {
        Configuration cfg = context.getConfiguration();

        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues();
        // Run flags
        final boolean runFastQCOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false);
        final boolean runFastQC = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC, true);
        final boolean runAlignmentOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false);

        final boolean useExistingPairedBams = cfgValues.getBoolean(COConstants.FLAG_USE_EXISTING_PAIRED_BAMS, false);
        // Usage flags
        final boolean useCombinedAlignAndSampe = cfgValues.getBoolean(COConstants.FLAG_USE_COMBINED_ALIGN_AND_SAMPE, false);
        final boolean runSlimWorkflow = cfgValues.getBoolean(COConstants.FLAG_RUN_SLIM_WORKFLOW, false);

        BamFileGroup sortedBamFiles = new BamFileGroup();

        if (useExistingPairedBams) {
            //Start from the paired bams instead of the lane files.
            sortedBamFiles = runtimeService.getPairedBamFilesForDataSet(context, sample);
        } else {

            //Create bam files out of the lane files
            List<LaneFileGroup> rawSequenceGroups = loadLaneFilesForSample(context, sample);
            if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
                return sortedBamFiles;
            for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                if (runFastQC && !runAlignmentOnly)
                    rawSequenceGroup.calcFastqcForAll();
                if (runFastQCOnly)
                    continue;


                BamFile bamFile = null;

                if (useCombinedAlignAndSampe) { //I.e. bwa mem
                    if (runSlimWorkflow) {
                        bamFile = rawSequenceGroup.alignAndPairSlim();
                    } else {
                        bamFile = rawSequenceGroup.alignAndPair();
                    }
                } else { //I.e. bwa align
                    rawSequenceGroup.alignAll();
                    if (runSlimWorkflow) {
                        bamFile = rawSequenceGroup.getAllAlignedFiles().pairAndSortSlim();
                    } else {
                        bamFile = rawSequenceGroup.getAllAlignedFiles().pairAndSort();
                    }
                }

                bamFile.setAsTemporaryFile();  // Bam files created with sai files are only temporary.
                sortedBamFiles.addFile(bamFile);
            }

        }
        return sortedBamFiles;
    }

    private BamFile mergeAndRemoveDuplicatesFat(ExecutionContext context, Sample sample, BamFileGroup sortedBamFiles) {
        Configuration cfg = context.getConfiguration();
        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues();
        final boolean runCollectBamFileMetrics = cfgValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false);

        //To avoid problems with qcsummary the step is done manually.
        sortedBamFiles.runDefaultOperations();
        for (BamFile sortedBam : sortedBamFiles.getFilesInGroup())
            sortedBam.createQCSummaryFile();

        BamFile mergedBam = sortedBamFiles.mergeAndRemoveDuplicates();

        if (runCollectBamFileMetrics) mergedBam.collectMetrics();
        mergedBam.runDefaultOperations();
        mergedBam.calcCoverage();
        mergedBam.createQCSummaryFile();

        return mergedBam;
    }

    private boolean valueIsEmpty(ExecutionContext context, Object value, String variableName) {
        if (value == null || value.toString() == "") {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("Expected value to be set: ${variableName}"))
            return true
        }
        return false
    }

    private boolean fileIsAccessible(ExecutionContext context, File file, String variableName) {
        if (valueIsEmpty(context, file, variableName) || !FileSystemAccessProvider.getInstance().checkFile(file, false, context)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("File '${file}' not accessible: ${variableName}"))
            return false
        }
        return true
    }

    private boolean directoryIsAccessible(ExecutionContext context, File directory, String variableName) {
        if (valueIsEmpty(context, directory, variableName) || !FileSystemAccessProvider.getInstance().checkDirectory(directory, context, false)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("Directory '${directory}' not accessible: ${variableName}"))
            return false
        }
        return true
    }

    private boolean checkConfiguration(ExecutionContext context) {
        AlignmentAndQCConfig config = new AlignmentAndQCConfig(context)

        boolean returnValue
        returnValue =
                !valueIsEmpty(context, config.getIndexPrefix(), AlignmentAndQCConfig.CVALUE_INDEX_PREFIX) &&
                directoryIsAccessible(context, new File(config.getIndexPrefix()).getParentFile(), AlignmentAndQCConfig.CVALUE_INDEX_PREFIX)
        returnValue =
                fileIsAccessible(context, config.getChromosomeSizesFile(), AlignmentAndQCConfig.CVALUE_CHROMOSOME_SIZES_FILE)
        if (config.getRunExomeAnalysis()) {
            returnValue =
                    fileIsAccessible(context, config.getTargetRegionsFile(), AlignmentAndQCConfig.CVALUE_TARGET_REGIONS_FILE) &&
                    !valueIsEmpty(context, config.getTargetSize(), AlignmentAndQCConfig.CVALUE_TARGET_SIZE)
        }
        return returnValue
    }

    private boolean checkSamples(ExecutionContext context) {
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (samples.size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("No samples found for PID ${context.getDataSet()}!"))
            return false
        } else {
            logger.postAlwaysInfo("Found " + samples.size() + " samples for dataset " + context.getDataSet().getId());
            return true
        }
    }

    protected boolean checkLaneFiles(ExecutionContext context) {
        boolean returnValue = true
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService();
        List<Sample> samples = runtimeService.getSamplesForContext(context);
        final boolean useExistingPairedBams = context.getConfiguration().getConfigurationValues().getBoolean(COConstants.FLAG_USE_EXISTING_PAIRED_BAMS, false);
        if (!useExistingPairedBams) {
            int cnt = 0;
            for (Sample sample : samples) {
                List<LaneFileGroup> laneFileGroups = loadLaneFilesForSample(context, sample);
                for (LaneFileGroup lfg : laneFileGroups) {
                    cnt += lfg.getFilesInGroup().size();
                }
                logger.postAlwaysInfo("Processed sample " + sample.getName() + " and found " + laneFileGroups.size() + " groups of lane files.");
            }
            if (cnt <= 0) {
                context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                        expand("No lane files found for PID ${context.getDataSet()}!"))
                returnValue = false
            }
        }
        return returnValue;
    }

    protected boolean checkSingleBam(ExecutionContext context) {

        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        if (!aqcfg.getSingleBamParameter()) return true

        boolean returnValue = true
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (samples.size() > 1) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("A bam parameter for single bam was set, but there is more than one sample available."));
            returnValue &= false;
        }

        def accessProvider = FileSystemAccessProvider.getInstance()
        def bamFile = new File(aqcfg.getSingleBamParameter())
        if (!accessProvider.fileExists(bamFile) || !accessProvider.isReadable(bamFile)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("A bam parameter for single bam was set, but the bam file is not readable: '${bamFile}'"));
            returnValue &= false;
        }

        return returnValue
    }

    public boolean checkFingerprintingSitesFile(ExecutionContext context) {
        def aqcfg = new AlignmentAndQCConfig(context)
        def accessProvider = FileSystemAccessProvider.getInstance()
        boolean result = true
        if (aqcfg.runFingerprinting) {
            if (!accessProvider.fileExists(aqcfg.fingerprintingSitesFile)
                    || !accessProvider.isReadable(aqcfg.fingerprintingSitesFile)) {
                context.addErrorEntry(ExecutionContextError.
                        EXECUTION_SETUP_INVALID.expand("Fingerprinting reference sites file not readable: '${aqcfg.fingerprintingSitesFile}'"))
                result = false
            }
        }
        return result
    }

    @Override
    public boolean checkExecutability(ExecutionContext context) {
        boolean result = super.checkExecutability(context)
        result &= checkConfiguration(context)
        result &= checkSamples(context)
        result &= checkLaneFiles(context)
        result &= checkSingleBam(context)
        result &= checkFingerprintingSitesFile(context)
        return result
    }


    @Override
    public boolean createTestdata(ExecutionContext context) {
        boolean allOk = true;
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService();

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
