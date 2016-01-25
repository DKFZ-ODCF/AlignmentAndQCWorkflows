package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.common.BasicCOProjectsRuntimeService;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.StringConstants;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.core.ExecutionContextError

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES;

/**
 * @author michael
 */
public class QCPipelinePancancerSlim extends QCPipeline {

    @Override
    public boolean execute(ExecutionContext context) {
        Configuration cfg = context.getConfiguration();
        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues();
        cfgValues.put(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, "false", "boolean"); //Disable sample extraction from output for alignment workflows.

        // Run flags
        final boolean runAlignmentOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false);
        final boolean runCoveragePlots = cfgValues.getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true);
        final boolean runExomeAnalysis = cfgValues.getBoolean(COConstants.FLAG_RUN_EXOME_ANALYSIS);
        final boolean runCollectBamFileMetrics = cfgValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false);
        final boolean runCoveragePlotsOnly = cfgValues.getBoolean("runCoveragePlotsOnly", false);

        final String overrideFastqFiles = cfgValues.getString("overrideFastqFiles", "");
        final String[] overrideMergedBamFiles = cfgValues.getString("overrideBamFiles", "").split(StringConstants.SPLIT_SEMICOLON);
        final String[] overrideSampleNames = cfgValues.getString("overrideSampleNames", "").split(StringConstants.SPLIT_SEMICOLON);
        String controlBamName = overrideMergedBamFiles[0];
        String[] tumorBamNames = overrideMergedBamFiles[1 .. -1];

        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getProject().getRuntimeService();

        List<Sample> samples
        if(overrideSampleNames) {
            samples = overrideSampleNames.collect { String sname -> new Sample(context, sname) }
        } else {
            samples = runtimeService.getSamplesForContext(context);
        }
        if (!samples)
            return false;

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = [:]

        for (Sample sample in samples) {
            BamFileGroup sortedBamFiles = createSortedBams(context, runtimeService, sample);

            if (!sortedBamFiles.getFilesInGroup()) continue;

            if (runAlignmentOnly) continue;

            BamFile mergedBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample);
            if (runCollectBamFileMetrics) mergedBam.collectMetrics();

            if (runExomeAnalysis)
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

        if (runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot();
        }

        return true;
    }

    private BamFileGroup createSortedBams(ExecutionContext context, BasicCOProjectsRuntimeService runtimeService, Sample sample) {
        Configuration cfg = context.getConfiguration();
        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues();
        // Run flags
        final boolean runFastQCOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false);
        final boolean runFastQC = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC, true);
        final boolean runAlignmentOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false);

        BamFileGroup sortedBamFiles = new BamFileGroup();

        //Create bam files out of the lane files
        List<LaneFileGroup> rawSequenceGroups = loadLaneFilesForSample(context, sample);
        if (rawSequenceGroups == null || rawSequenceGroups.size() == 0)
            return sortedBamFiles;
        for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
            if (runFastQC && !runAlignmentOnly)
                rawSequenceGroup.calcFastqcForAll();
            if (runFastQCOnly)
                continue;

            BamFile bamFile = rawSequenceGroup.alignAndPairSlim();

            bamFile.setAsTemporaryFile();  // Bam files created with sai files are only temporary.
            sortedBamFiles.addFile(bamFile);
        }

        return sortedBamFiles;
    }
}
