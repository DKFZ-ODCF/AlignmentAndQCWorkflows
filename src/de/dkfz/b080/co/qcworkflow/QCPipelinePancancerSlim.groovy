package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.common.AlignmentAndQCConfig
import de.dkfz.b080.co.common.BasicCOProjectsRuntimeService
import de.dkfz.b080.co.files.*
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES

class QCPipelinePancancerSlim extends QCPipeline {

    @Override
    boolean execute(ExecutionContext context) {
        // Disable sample extraction from output for alignment workflows.
        context.getConfiguration().getConfigurationValues()
                .put(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, "false", "boolean")

        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        BasicCOProjectsRuntimeService runtimeService =
                (BasicCOProjectsRuntimeService) context.runtimeService

        List<Sample> samples
        if(!aqcfg.overrideSampleNames.empty) {
            samples = aqcfg.overrideSampleNames.collect { String sname -> new Sample(context, sname) }
        } else {
            samples = runtimeService.metadataAccessor.getSamples(context)
        }
        if (!samples)
            return false

        BamFileGroup mergedBamFiles = new BamFileGroup()
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = [:]

        for (Sample sample in samples) {
            BamFileGroup sortedBamFiles = createSortedBams(context, runtimeService, sample)

            if (!sortedBamFiles.getFilesInGroup()) continue

            if (aqcfg.runAlignmentOnly) continue

            BamFile mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample)

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

        if (!mergedBamFiles.getFilesInGroup()) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."))
            return false
        }

        if (aqcfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR))
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot()
        }

        return true
    }

    private BamFileGroup createSortedBams(ExecutionContext context, BasicCOProjectsRuntimeService runtimeService, Sample sample) {
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        BamFileGroup sortedBamFiles = new BamFileGroup()

        // Create bam files out of the lane files
        List<LaneFileGroup> rawSequenceGroups = loadLaneFilesForSample(context, sample)
        if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
            return sortedBamFiles

        for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
            if (aqcfg.runFastqc && !aqcfg.runAlignmentOnly)
                rawSequenceGroup.calcFastqcForAll()
            if (aqcfg.runFastqcOnly)
                continue

            BamFile bamFile = rawSequenceGroup.alignAndPairSlim()

            bamFile.setAsTemporaryFile()  // Bam files created with sai files are only temporary.
            sortedBamFiles.addFile(bamFile)
        }

        return sortedBamFiles
    }
}
