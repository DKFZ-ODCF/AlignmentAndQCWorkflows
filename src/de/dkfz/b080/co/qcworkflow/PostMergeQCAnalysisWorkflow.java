package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.common.AlignmentAndQCConfig;
import de.dkfz.b080.co.common.COProjectsRuntimeService;
import de.dkfz.b080.co.files.BamFile;
import de.dkfz.b080.co.files.BamFileGroup;
import de.dkfz.b080.co.files.CoverageTextFileGroup;
import de.dkfz.b080.co.files.Sample;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.core.ExecutionContextError;
import de.dkfz.roddy.core.Workflow;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * A short workflow which only does post merge alignment quality control
 */
public class PostMergeQCAnalysisWorkflow extends Workflow {

    public PostMergeQCAnalysisWorkflow() {}

    @Override
    public boolean execute(ExecutionContext context) {
        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context);
        cfg.setExtractSamplesFromOutputFiles(true); // Enable sample extraction from output for this workflow.
        cfg.setUseOnlyExistingTargetBam(true);  // Enable usage of existing merged bams for this workflow.

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService();
        List<Sample> samples = runtimeService.metadataAccessor.getSamples(context);
        if (samples.size() == 0)
            return false;

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>();

        for (Sample sample : samples) {
            BamFile mergedBam = new BamFile(runtimeService.metadataAccessor.
                    getMergedBamFileFromFilesystem(context, null, sample));

            mergedBam.performPostMergeQCAnalysis();
            if (cfg.getRunExomeAnalysis()) {
                mergedBam.extractTargetsCalculateCoverage();
            }

            Sample.SampleType sampleType = sample.getType();
            if (!coverageTextFilesBySample.containsKey(sampleType))
                coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup());
            coverageTextFilesBySample.get(sampleType).addFile(mergedBam.getReadBinsCoverageTextFile());

            mergedBamFiles.addFile(mergedBam);
        }

        if (mergedBamFiles.getFilesInGroup().size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."));
            return false;
        }

        if (cfg.getRunCoveragePlots() && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot();
        }

        return true;
    }
}
