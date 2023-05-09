package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.files.*
import de.dkfz.b080.co.common.*
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.tools.LoggerWrapper
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues
import de.dkfz.roddy.core.*
import groovy.transform.CompileStatic

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES

/**
 * @author michael
 */
@CompileStatic
class QCPipeline extends Workflow {

    private static LoggerWrapper logger = LoggerWrapper.getLogger(QCPipeline.class.getName())

    QCPipeline() {}

    @Override
    boolean execute(ExecutionContext context) {
        // Disable sample extraction from output for alignment workflows.
        context.getConfiguration().configurationValues.
                put(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, "false", "boolean")

        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        COProjectsRuntimeService runtimeService =
                (COProjectsRuntimeService) context.runtimeService

        List<Sample> samples = runtimeService.metadataAccessor.getSamples(context)
        if (samples.size() == 0)
            return false

        BamFileGroup mergedBamFiles = new BamFileGroup()
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>()

        for (Sample sample : samples) {
            BamFileGroup sortedBamFiles = createSortedBams(context, runtimeService, sample)

            if (sortedBamFiles.getFilesInGroup().size() == 0) continue

            if (aqcfg.runFastqcOnly || aqcfg.runAlignmentOnly) continue

            BamFile mergedBam
            mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample)
            
            if (aqcfg.runCollectBamFileMetrics)
                mergedBam.collectMetrics()
            
            if (aqcfg.runExomeAnalysis)
                mergedBam.extractTargetsCalculateCoverage()

            Sample.SampleType sampleType = sample.getType()
            if (!coverageTextFilesBySample.containsKey(sampleType))
                coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup())
            coverageTextFilesBySample.get(sampleType).addFile(mergedBam.calcReadBinsCoverage())

            mergedBamFiles.addFile(mergedBam)
        }

        if (mergedBamFiles.getFilesInGroup().size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."))
            return false
        }

        if (aqcfg.runFastqcOnly)
            return true

        if (aqcfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR))
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot()
        }

        return true
    }

    private BamFileGroup createSortedBams(ExecutionContext context, COProjectsRuntimeService runtimeService, Sample sample) {
        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context)
        BamFileGroup sortedBamFiles = new BamFileGroup()

        if (cfg.useOnlyExistingPairedBams) {
            //Start from the paired bams instead of the lane files.
            sortedBamFiles = runtimeService.getPairedBamFilesForDataSet(context, sample)
        } else {

            //Create bam files out of the lane files
            List<LaneFileGroup> rawSequenceGroups = ((COProjectsRuntimeService)context.getRuntimeService()).loadLaneFilesForSample(context, sample)
            if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
                return sortedBamFiles
            for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                if (cfg.runFastqc && !cfg.runAlignmentOnly)
                    rawSequenceGroup.calcFastqcForAll()
                if (cfg.runFastqcOnly)
                    continue
                
                BamFile bamFile = rawSequenceGroup.alignAndPairSlim()

                bamFile.setAsTemporaryFile();  // Bam files created with sai files are only temporary.
                sortedBamFiles.addFile(bamFile)
            }

        }
        return sortedBamFiles
    }

    private BamFile mergeAndRemoveDuplicatesFat(ExecutionContext context, Sample sample, BamFileGroup sortedBamFiles) {
        Configuration cfg = context.getConfiguration()
        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues()
        final boolean runCollectBamFileMetrics = cfgValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false)

        //To avoid problems with qcsummary the step is done manually.
        sortedBamFiles.runDefaultOperations()
        for (BamFile sortedBam : sortedBamFiles.filesInGroup)
            sortedBam.createQCSummaryFile()

        BamFile mergedBam = sortedBamFiles.mergeAndRemoveDuplicates()

        if (runCollectBamFileMetrics) mergedBam.collectMetrics()
        mergedBam.runDefaultOperations()
        mergedBam.calcCoverage()
        mergedBam.createQCSummaryFile()

        return mergedBam
    }

    private boolean valueIsEmpty(ExecutionContext context, Object value, String variableName) {
        if (value == null || value.toString() == "") {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("Expected value to be set: ${variableName}"))
            return true
        }
        return false
    }

    private boolean fileIsAccessible(ExecutionContext context, File file, String variableName) {
        if (valueIsEmpty(context, file, variableName) ||
                !FileSystemAccessProvider.instance.checkFile(file, false, context)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("File '${file}' not accessible: ${variableName}"))
            return false
        }
        return true
    }

    private boolean directoryIsAccessible(ExecutionContext context, File directory, String variableName) {
        if (valueIsEmpty(context, directory, variableName) ||
                !FileSystemAccessProvider.instance.checkDirectory(directory, context, false)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("Directory '${directory}' not accessible: ${variableName}"))
            return false
        }
        return true
    }

    private boolean checkConfiguration(ExecutionContext context) {
        AlignmentAndQCConfig config = new AlignmentAndQCConfig(context)

        boolean returnValue
        returnValue =
                !valueIsEmpty(context, config.indexPrefix,
                        AlignmentAndQCConfig.CVALUE_INDEX_PREFIX) &&
                directoryIsAccessible(context, new File(config.indexPrefix).parentFile,
                        AlignmentAndQCConfig.CVALUE_INDEX_PREFIX)
        returnValue &=
                fileIsAccessible(context, config.chromosomeSizesFile,
                        AlignmentAndQCConfig.CVALUE_CHROMOSOME_SIZES_FILE)
        if (config.runExomeAnalysis) {
            returnValue &=
                    fileIsAccessible(context, config.targetRegionsFile,
                            AlignmentAndQCConfig.CVALUE_TARGET_REGIONS_FILE) &&
                    !valueIsEmpty(context, config.targetSize,
                            AlignmentAndQCConfig.CVALUE_TARGET_SIZE)
        }
        return returnValue
    }

    private boolean checkSamples(ExecutionContext context) {
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.runtimeService
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (samples.size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("No samples found for PID ${context.dataSet}!"))
            return false
        } else {
            logger.postAlwaysInfo("Found " + samples.size() + " samples for dataset " + context.dataSet.id)
            return true
        }
    }

    protected boolean checkLaneFiles(ExecutionContext context) {
        boolean returnValue = true

        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context)

        BasicCOProjectsRuntimeService runtimeService =
                (BasicCOProjectsRuntimeService) context.runtimeService
        List<Sample> samples = runtimeService.getSamplesForContext(context)

        if (cfg.useOnlyExistingTargetBam && cfg.fastqFileListIsSet) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                    expand("Both 'fastq_list' and 'useOnlyExistingTargetBam' are set. Set only one of them!"))
            returnValue = false
        }

        if (cfg.useOnlyExistingTargetBam && cfg.useExistingLaneBams) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("Both 'useExistingLaneBams' and 'useOnlyExistingTargetBam' are set. Set only one of them!"))
            returnValue = false
        }

        if (!cfg.useExistingLaneBams) {
            int cnt = 0
            for (Sample sample : samples) {
                List<LaneFileGroup> laneFileGroups = ((COProjectsRuntimeService) runtimeService).loadLaneFilesForSample(context, sample)
                for (LaneFileGroup lfg : laneFileGroups) {
                    cnt += lfg.filesInGroup.size()
                }
                logger.postAlwaysInfo("Processed sample " + sample.name + " and found " + laneFileGroups.size() + " groups of lane files.")
            }
            if (cnt <= 0) {
                context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                        expand("No lane files found for PID ${context.dataSet}!"))
                returnValue = false
            }
        }

        return returnValue
    }

    protected boolean checkSingleBam(ExecutionContext context) {

        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        if (!aqcfg.singleBamParameter) return true

        boolean returnValue = true

        if (aqcfg.useOnlyExistingTargetBam) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("Both 'bam' and 'useOnlyExistingTargetBam' are set. Set only one of them!"))
            returnValue &= false
        }

        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (samples.size() > 1) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("A 'bam' parameter for a single BAM, but there is more than one sample available."))
            returnValue &= false
        }

        def accessProvider = FileSystemAccessProvider.getInstance()
        def bamFile = new File(aqcfg.singleBamParameter)
        if (!accessProvider.fileExists(bamFile) || !accessProvider.isReadable(bamFile)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("A 'bam' parameter was set, but the BAM file is not readable: '${bamFile}'"))
            returnValue &= false
        }

        return returnValue
    }

    boolean checkFingerprintingSitesFile(ExecutionContext context) {
        def aqcfg = new AlignmentAndQCConfig(context)
        def accessProvider = FileSystemAccessProvider.instance
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
    boolean checkExecutability(ExecutionContext context) {
        boolean result = super.checkExecutability(context)
        result &= checkConfiguration(context)
        result &= checkSamples(context)
        result &= checkLaneFiles(context)
        result &= checkSingleBam(context)
        result &= checkFingerprintingSitesFile(context)
        return result
    }


    boolean createTestdata(ExecutionContext context) {
        boolean allOk = true
        COProjectsRuntimeService runtimeService =
                (COProjectsRuntimeService) context.runtimeService

        List<Sample> samples = runtimeService.getSamplesForContext(context)
        for (Sample sample : samples) {
            List<LaneFile> files = new LinkedList<LaneFile>()
            LaneFileGroup allLaneFiles = new LaneFileGroup(context, "allLaneFiles", "noSpecificRun", sample, files)

            List<LaneFileGroup> rawSequenceGroups = runtimeService.getLanesForSample(context, sample)
            for (LaneFileGroup lfg : rawSequenceGroups) {
                for (LaneFile lf : lfg.filesInGroup) {
                    allLaneFiles.addFile(lf)
                }
            }
            allLaneFiles.createTestDataForLaneFiles()
        }
        return allOk
    }
}
