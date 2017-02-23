package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.files.*;
import de.dkfz.b080.co.common.*
import de.dkfz.roddy.execution.io.ExecutionService
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
    public boolean execute(ExecutionContext context) { ;
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        aqcfg.extractSamplesFromOutputFiles = false

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForContext(context);
        if (samples.size() == 0)
            return false;

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>();

        for (Sample sample : samples) {
            BamFileGroup sortedBamFiles = createSortedBams(context, runtimeService, sample);

            if (sortedBamFiles.getFilesInGroup().size() == 0) continue;

            if (aqcfg.runFastQCOnly || aqcfg.runAlignmentOnly) continue;

            BamFile mergedBam;
            if (aqcfg.runSlimWorkflow) {
                mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample);
                if (aqcfg.runCollectBamFileMetrics) mergedBam.collectMetrics();
            } else {
                mergedBam = mergeAndRemoveDuplicatesFat(context, sample, sortedBamFiles);
            }

            if (aqcfg.runExomeAnalysis) {
                if (!aqcfg.runSlimWorkflow) mergedBam.rawBamCoverage();
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

        if (aqcfg.runFastQCOnly)
            return true;

        if (aqcfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot();
        }

        return true;
    }

    private BamFileGroup createSortedBams(ExecutionContext context, COProjectsRuntimeService runtimeService, Sample sample) {
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        BamFileGroup sortedBamFiles = new BamFileGroup();

        if (aqcfg.useExistingPairedBams) {
            //Start from the paired bams instead of the lane files.
            sortedBamFiles = runtimeService.getPairedBamFilesForDataSet(context, sample);
        } else {

            //Create bam files out of the lane files
            List<LaneFileGroup> rawSequenceGroups = ((COProjectsRuntimeService)context.getRuntimeService()).loadLaneFilesForSample(context, sample);
            if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
                return sortedBamFiles;
            for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                if (aqcfg.runFastQC && !aqcfg.runAlignmentOnly)
                    rawSequenceGroup.calcFastqcForAll();
                if (aqcfg.runFastQCOnly)
                    continue;


                BamFile bamFile = null;

                if (aqcfg.useCombinedAlignAndSampe) { //I.e. bwa mem
                    if (aqcfg.runSlimWorkflow) {
                        bamFile = rawSequenceGroup.alignAndPairSlim();
                    } else {
                        bamFile = rawSequenceGroup.alignAndPair();
                    }
                } else { //I.e. bwa align
                    rawSequenceGroup.alignAll();
                    if (aqcfg.runSlimWorkflow) {
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
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)

        //To avoid problems with qcsummary the step is done manually.
        sortedBamFiles.runDefaultOperations();
        for (BamFile sortedBam : sortedBamFiles.getFilesInGroup())
            sortedBam.createQCSummaryFile();

        BamFile mergedBam = sortedBamFiles.mergeAndRemoveDuplicates();

        if (aqcfg.runCollectBamFileMetrics) mergedBam.collectMetrics();
        mergedBam.runDefaultOperations();
        mergedBam.calcCoverage();
        mergedBam.createQCSummaryFile();

        return mergedBam;
    }

    boolean valueIsEmpty(ExecutionContext context, Object value, String variableName) {
        if (value == null || value.toString() == "") {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("Expected value to be set: ${variableName}"))
            return true
        }
        return false
    }

    boolean fileIsAccessible(ExecutionContext context, File file, String variableName) {
        if (valueIsEmpty(context, file, variableName) || !FileSystemAccessProvider.getInstance().checkFile(file, false, context)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("File '${file}' not accessible: ${variableName}"))
            return false
        }
        return true
    }

    boolean directoryIsAccessible(ExecutionContext context, File directory, String variableName) {
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
        returnValue &=
                fileIsAccessible(context, config.getChromosomeSizesFile(), AlignmentAndQCConfig.CVALUE_CHROMOSOME_SIZES_FILE)
        if (config.getRunExomeAnalysis()) {
            returnValue &=
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
                List<LaneFileGroup> laneFileGroups = ((COProjectsRuntimeService)runtimeService).loadLaneFilesForSample(context, sample);
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

    /** Check that the fingerprinting sites file is accessible, unless it is in the plugin.
     *
     * @param context
     * @return
     */
    public boolean checkFingerprintingSitesFile(ExecutionContext context) {
        def aqcfg = new AlignmentAndQCConfig(context)
        boolean result = true
        if (aqcfg.runFingerprinting) {
            if (!aqcfg.fingerprintingSitesFile.absolutePath.contains("\$${ExecutionService.RODDY_CVALUE_DIRECTORY_EXECUTION}")) {
                result = fileIsAccessible(context, aqcfg.fingerprintingSitesFile, AlignmentAndQCConfig.CVALUE_FINGERPRINTING_SITES_FILE)
            }
        }
        return result
    }

    /** Check that the clipping index file is accessible, unless it is in the plugin.
     *
     * @param context
     * @return
     */
    public boolean checkClipIndexFile(ExecutionContext context) {
        def aqcfg = new AlignmentAndQCConfig(context)
        Boolean result = true
        if (aqcfg.useAdapterTrimming) {
            if (!aqcfg.clipIndex.absolutePath.contains("\$${ExecutionService.RODDY_CVALUE_DIRECTORY_EXECUTION}")) {
                result = fileIsAccessible(context, aqcfg.clipIndex, AlignmentAndQCConfig.CVALUE_CLIP_INDEX)
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
        result &= checkClipIndexFile(context)
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
