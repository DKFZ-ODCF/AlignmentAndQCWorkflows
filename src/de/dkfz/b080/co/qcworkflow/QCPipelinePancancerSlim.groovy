package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.common.AlignmentAndQCConfig;
import de.dkfz.b080.co.common.BasicCOProjectsRuntimeService;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.core.ExecutionContextError

/**
 * @author michael
 */
public class QCPipelinePancancerSlim extends QCPipeline {

    @Override
    public boolean execute(ExecutionContext context) {
        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context)
        cfg.sampleExtractionFromOutputFiles = false

        String controlBamName = overrideMergedBamFiles[0];
        String[] tumorBamNames = overrideMergedBamFiles[1 .. -1];

        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService();

        List<Sample> samples
        if(cfg.overrideSampleNames) {
            samples = cfg.overrideSampleNames.collect { String sname -> new Sample(context, sname) }
        } else {
            samples = runtimeService.getSamplesForContext(context);
        }
        if (!samples)
            return false;

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = [:]

        for (Sample sample in samples) {
            BamFileGroup sortedBamFiles = createSortedBams(cfg, runtimeService, sample);

            if (!sortedBamFiles.getFilesInGroup()) continue;

            if (cfg.runAlignmentOnly) continue;

            BamFile mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample);
            if (cfg.runCollectBamFileMetrics) mergedBam.collectMetrics();

            if (cfg.runExomeAnalysis)
                BamFile targetOnlyBamFile = mergedBam.extractTargetsCalculateCoverage();

            Sample.SampleType sampleType = sample.getType();
            if (!coverageTextFilesBySample.containsKey(sampleType))
                coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup());
            coverageTextFilesBySample.get(sampleType).addFile(mergedBam.calcReadBinsCoverage());

            mergedBamFiles.addFile(mergedBam);
        }

        if (!mergedBamFiles.getFilesInGroup()) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."));
            return false;
        }

        if (cfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot();
        }

        return true;
    }

    private BamFileGroup createSortedBams(ExecutionContext context, BasicCOProjectsRuntimeService runtimeService, Sample sample) {
        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context)
        BamFileGroup sortedBamFiles = new BamFileGroup();

        //Create bam files out of the lane files
        List<LaneFileGroup> rawSequenceGroups = loadLaneFilesForSample(cfg.context, sample);
        if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
            return sortedBamFiles;
        for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
            if (cfg.runFastqQC && !cfg.runAlignmentOnly)
                rawSequenceGroup.calcFastqcForAll();
            if (cfg.runFastQCOnly)
                continue;

            BamFile bamFile = rawSequenceGroup.alignAndPairSlim();
            bamFile.setAsTemporaryFile();  // Bam files created with sai files are only temporary.
            sortedBamFiles.addFile(bamFile);
        }

        return sortedBamFiles;
    }
}
