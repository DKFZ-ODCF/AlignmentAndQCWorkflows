/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.common.AlignmentAndQCConfig
import de.dkfz.b080.co.common.BasicCOProjectsRuntimeService
import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.*
import de.dkfz.b080.co.methods.ACEseq
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.core.Workflow
import de.dkfz.roddy.execution.io.ExecutionService
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.tools.LoggerWrapper

/**
 * @author michael
 */
@groovy.transform.CompileStatic
class QCPipeline extends Workflow {

    private static LoggerWrapper logger = LoggerWrapper.getLogger(QCPipeline.class.getName())

    QCPipeline() {}

    @Override
    boolean execute(ExecutionContext context) {
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        aqcfg.extractSamplesFromOutputFiles = false

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()

        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (samples.size() == 0)
            return false

        BamFileGroup mergedBamFiles = new BamFileGroup()
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>()

        for (Sample sample : samples) {
            BamFileGroup sortedBamFiles = createLaneBams(aqcfg, runtimeService, sample)

            if (sortedBamFiles.getFilesInGroup().size() == 0) continue

            if (aqcfg.runFastQCOnly || aqcfg.runAlignmentOnly) continue

            BamFile mergedBam
            mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample)
            if (aqcfg.runCollectBamFileMetrics) mergedBam.collectMetrics()

            if (aqcfg.runExomeAnalysis) {
                mergedBam.rawBamCoverage()
                mergedBam.extractTargetsCalculateCoverage()
            }


            if (!coverageTextFilesBySample.containsKey(sample.sampleType))
                coverageTextFilesBySample.put(sample.sampleType, new CoverageTextFileGroup())
            coverageTextFilesBySample.get(sample.sampleType).addFile(mergedBam.calcReadBinsCoverage())

            mergedBamFiles.addFile(mergedBam)

            // The ACEseq QC could also be done per lane/run, but for non X10 data there is not the required >30x coverage (Kortine).
            if (aqcfg.runACEseqQC) {
                ACEseq.aceSeqQc(mergedBam.readBinsCoverageTextFile, sample)
            }
        }

        if (mergedBamFiles.getFilesInGroup().size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."))
            return false
        }

        if (aqcfg.runFastQCOnly)
            return true

        if (aqcfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR))
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot()
        }

        return true
    }

    private BamFileGroup createLaneBams(AlignmentAndQCConfig aqcfg, COProjectsRuntimeService runtimeService, Sample sample) {
        BamFileGroup sortedBamFiles = new BamFileGroup()

        if (aqcfg.useOnlyExistingLaneBams) {
            // Start from the paired BAMs instead of the FASTQ files.
            sortedBamFiles = runtimeService.getPairedBamFilesForDataSet(aqcfg.context, sample)

        } else {
            // Create bam files out of the FASTQ files.
            List<LaneFileGroup> rawSequenceGroups = runtimeService.loadLaneFilesForSample(aqcfg.context, sample)
            if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
                return sortedBamFiles
            for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                if (aqcfg.runFastQC && !aqcfg.runAlignmentOnly)
                    rawSequenceGroup.calcFastqcForAll()
                if (aqcfg.runFastQCOnly)
                    continue

                BamFile bamFile = null

                if (aqcfg.useCombinedAlignAndSampe) { //I.e. bwa mem
                    bamFile = rawSequenceGroup.alignAndPairSlim()
                } else { //I.e. bwa align
                    rawSequenceGroup.alignAll()
                    bamFile = rawSequenceGroup.getAllAlignedFiles().pairAndSortSlim()
                }

                // @Michael: The comment suggests that this should only be called in the else branch above!?
                // @Michael: Why should BAM files created with sai files be temporary? Because they are lane BAMs?
                bamFile.setAsTemporaryFile()  // Bam files created with sai files are only temporary.
                sortedBamFiles.addFile(bamFile)
            }

        }
        return sortedBamFiles
    }

    private static ExecutionContextError getInvalidError (String message) {
        return ExecutionContextError.EXECUTION_SETUP_INVALID.expand(message)
    }

    private boolean checkConfiguration(ExecutionContext context) {
        AlignmentAndQCConfig config = new AlignmentAndQCConfig(context)
        boolean returnValue
        returnValue =
                !context.valueIsEmpty(config.getIndexPrefix(), AlignmentAndQCConfig.CVALUE_INDEX_PREFIX) &&
                context.directoryIsAccessible(new File(config.getIndexPrefix()).getParentFile(), AlignmentAndQCConfig.CVALUE_INDEX_PREFIX)
        returnValue &=
                context.fileIsAccessible(config.getChromosomeSizesFile(), AlignmentAndQCConfig.CVALUE_CHROMOSOME_SIZES_FILE)
        if (config.runExomeAnalysis) {
            returnValue &=
                    context.fileIsAccessible(config.getTargetRegionsFile(), AlignmentAndQCConfig.CVALUE_TARGET_REGIONS_FILE) &&
                    !context.valueIsEmpty(config.getTargetSize(), AlignmentAndQCConfig.CVALUE_TARGET_SIZE)
        }
        if (config.useOnlyExistingTargetBam && config.fastqFileListIsSet) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                    expand("Both 'fastq_list' and 'useOnlyExistingTargetBam' are set. Set only one of them!"))
            returnValue = false
        }
        if (config.useOnlyExistingTargetBam && config.useOnlyExistingLaneBams) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("Both 'useExistingLaneBams' and 'useOnlyExistingTargetBam' are set. Set only one of them!"))
            returnValue = false
        }
        return returnValue
    }

    private boolean checkSamples(ExecutionContext context) {
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (samples.size() == 0) {
            context.addErrorEntry(getInvalidError("No samples found for PID ${context.getDataSet()}!"))
            return false
        } else {
            logger.postAlwaysInfo("Found " + samples.size() + " samples for dataset " + context.getDataSet().getId())
            return true
        }
    }

    private boolean checkFile(ExecutionContext context, File file) {
        if (!context.fileIsAccessible(file)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                    expand("Could not access '${file}'."))
            return false
        }
        return true
    }

    private <T extends BaseFile> boolean checkFiles(ExecutionContext context, List<T> files) {
        return files.collect { checkFile(context, it.path) }.inject { acc, i -> acc && i }
    }

    private Boolean checkLaneFileGroups(ExecutionContext context, List<LaneFileGroup> laneFileGroups) {
        boolean returnValue = true
        for (laneFileGroup in laneFileGroups) {
            if (new AlignmentAndQCConfig(context).useSingleEndProcessing)
                returnValue &= checkFiles(context, [laneFileGroup.filesInGroup[0]] as List)
            else
                returnValue &= checkFiles(context, laneFileGroup.filesInGroup)
        }
        return returnValue
    }

    protected boolean checkLaneFiles(ExecutionContext context) {
        boolean returnValue = true
        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context)
        COProjectsRuntimeService runtimeService = context.runtimeService as COProjectsRuntimeService
        if (!cfg.useOnlyExistingTargetBam) {
            for (Sample sample in runtimeService.getSamplesForContext(context)) {
                List<String> libraries = sample.getLibraries()
                if (libraries)
                    for (String library in libraries)
                        returnValue &= checkLaneFileGroups(context, runtimeService.loadLaneFilesForSampleAndLibrary(context, sample, library))
                else
                    returnValue &= checkLaneFileGroups(context, runtimeService.loadLaneFilesForSample(context, sample))
            }
        }
        return returnValue
    }

    protected boolean checkLaneFileCounts(ExecutionContext context) {
        boolean returnValue = true
        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context)
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        if (!cfg.useOnlyExistingLaneBams) {
            int cnt = 0
            for (Sample sample : samples) {
                List<LaneFileGroup> laneFileGroups = runtimeService.loadLaneFilesForSample(context, sample)
                for (LaneFileGroup lfg : laneFileGroups) {
                    cnt += lfg.getFilesInGroup().size()
                }
                logger.postAlwaysInfo("Processed sample " + sample.getName() + " and found " + laneFileGroups.size() + " groups of lane files.")
            }
            if (cnt <= 0) {
                context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                        expand("No lane files found for PID ${context.getDataSet()}!"))
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
        def bamFile = new File(aqcfg.getSingleBamParameter())
        if (!accessProvider.fileExists(bamFile) || !accessProvider.isReadable(bamFile)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("A 'bam' parameter was set, but the BAM file is not readable: '${bamFile}'"))
            returnValue &= false
        }

        return returnValue
    }

    /** Check that the fingerprinting sites file is accessible, unless it is in the plugin.
     *
     * @param context
     * @return
     */
    boolean checkFingerprintingSitesFile(ExecutionContext context) {
        def aqcfg = new AlignmentAndQCConfig(context)
        boolean result = true
        if (aqcfg.runFingerprinting) {
            if (!aqcfg.fingerprintingSitesFile.absolutePath.contains("\$${ExecutionService.RODDY_CVALUE_DIRECTORY_EXECUTION}")) {
                result = context.fileIsAccessible(aqcfg.fingerprintingSitesFile, AlignmentAndQCConfig.CVALUE_FINGERPRINTING_SITES_FILE)
            }
        }
        return result
    }

    /** If the ACEseq QC is turned on, check that the mappability, replication-timing and GC-content files are available.
     *
     * @param context
     * @return
     */
    boolean checkACEseqQCPrerequisites(ExecutionContext context) {
        def aqcfg = new AlignmentAndQCConfig(context)
        boolean result = true
        if (aqcfg.runACEseqQC) {
            // Kortine: The mappability file may have other than the same window size as the input data (i.e. WINDOW_SIZE),
            //          the replication timing and GC content files need to have the same window size as the input.
            if (aqcfg.mappabilityFile == null) {
                aqcfg.context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                        expand("The ACEseq QC steps require ${AlignmentAndQCConfig.CVALUE_MAPPABILITY_FILE} to be set."))
                result = false
            }
            result &= context.fileIsAccessible(aqcfg.mappabilityFile, AlignmentAndQCConfig.CVALUE_MAPPABILITY_FILE)

            if (aqcfg.replicationTimeFile == null) {
                aqcfg.context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                        expand("The ACEseq QC steps require ${AlignmentAndQCConfig.CVALUE_REPLICATION_TIME_FILE} to be set."))
                result = false
            }
            result &= context.fileIsAccessible(aqcfg.replicationTimeFile, AlignmentAndQCConfig.CVALUE_REPLICATION_TIME_FILE)

            if (aqcfg.gcContentFile == null) {
                aqcfg.context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                        expand("The ACEseq QC steps require ${AlignmentAndQCConfig.CVALUE_GC_CONTENT_FILE} to be set."))
                result = false
            }
            result &= context.fileIsAccessible(aqcfg.gcContentFile, AlignmentAndQCConfig.CVALUE_GC_CONTENT_FILE)
        }
        return result
    }

    /** Check that the clipping index file is accessible, unless it is in the plugin.
     *
     * @param context
     * @return
     */
    boolean checkClipIndexFile(ExecutionContext context) {
        def aqcfg = new AlignmentAndQCConfig(context)
        Boolean result = true
        if (aqcfg.useAdapterTrimming) {
            if (!aqcfg.clipIndex.absolutePath.contains("\$${ExecutionService.RODDY_CVALUE_DIRECTORY_EXECUTION}")) {
                result = context.fileIsAccessible(aqcfg.clipIndex, AlignmentAndQCConfig.CVALUE_CLIP_INDEX)
            }
        }
        return result
    }

    @Override
    boolean checkExecutability(ExecutionContext context) {
        boolean result = super.checkExecutability(context)
        result &= checkConfiguration(context)
        result &= checkLaneFileCounts(context)
        result &= checkLaneFiles(context)
        result &= checkSamples(context)
        result &= checkSingleBam(context)
        result &= checkFingerprintingSitesFile(context)
        result &= checkClipIndexFile(context)
        result &= checkACEseqQCPrerequisites(context)
        return result
    }

    boolean createTestdata(ExecutionContext context) {
        boolean allOk = true
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()

        List<Sample> samples = runtimeService.getSamplesForContext(context)
        for (Sample sample : samples) {
            List<LaneFile> files = new LinkedList<LaneFile>()
            LaneFileGroup allLaneFiles = new LaneFileGroup(context, "allLaneFiles", "noSpecificRun", sample, files)

            List<LaneFileGroup> rawSequenceGroups = runtimeService.getLanesForSample(context, sample)
            for (LaneFileGroup lfg : rawSequenceGroups) {
                for (LaneFile lf : lfg.getFilesInGroup()) {
                    allLaneFiles.addFile(lf)
                }
            }
            allLaneFiles.createTestDataForLaneFiles()
        }
        return allOk
    }
}
