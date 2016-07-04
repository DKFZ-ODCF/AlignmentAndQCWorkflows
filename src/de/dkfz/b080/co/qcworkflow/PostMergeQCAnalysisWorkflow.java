package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.common.COProjectsRuntimeService;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.core.ExecutionContextError;
import de.dkfz.roddy.core.Workflow;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES;
import static de.dkfz.b080.co.files.COConstants.FLAG_USE_EXISTING_MERGED_BAMS;

/**
 * A short workflow which only does post merge alignment quality control
 */
public class PostMergeQCAnalysisWorkflow extends Workflow {

    public PostMergeQCAnalysisWorkflow() {}

    @Override
    public boolean execute(ExecutionContext context) {
        Configuration cfg = context.getConfiguration();
        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues();
        cfgValues.put(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, "true", "boolean"); //Enable sample extraction from output for this workflow.
        cfgValues.put(FLAG_USE_EXISTING_MERGED_BAMS, "true", "boolean"); //Enable usage of existing merged bams for this workflow.
        final boolean runCoveragePlots = cfgValues.getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true);
        final boolean runExomeAnalysis = cfgValues.getBoolean(COConstants.FLAG_RUN_EXOME_ANALYSIS);

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getProject().getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForRun(context);
        if (samples.size() == 0)
            return false;

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>();

        for (Sample sample : samples) {
            BamFile mergedBam = runtimeService.getMergedBamFileForDataSetAndSample(context, sample);

            mergedBam.performPostMergeQCAnalysis();
            if (runExomeAnalysis) {
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

        if (runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (coverageTextFilesBySample.keySet().size() == 1) {
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot();
        }

        return true;
    }
}
