package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.common.BasicCOProjectsRuntimeService
import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.*
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues
import de.dkfz.roddy.core.DataSet
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError

import static de.dkfz.b080.co.files.COConstants.FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES

/**
 * @author michael
 */
@groovy.transform.CompileStatic
public class BisulfiteCoreWorkflow extends QCPipeline {

    @Override
    public boolean execute(ExecutionContext context) {
        Configuration cfg = context.getConfiguration();
        RecursiveOverridableMapContainerForConfigurationValues cfgValues = cfg.getConfigurationValues();
        cfgValues.put(FLAG_EXTRACT_SAMPLES_FROM_OUTPUT_FILES, "false", "boolean"); //Disable sample extraction from output for alignment workflows.

        // Run flags
        final boolean runFastQCOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false)
        final boolean runFastQC = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC, true)
        final boolean runAlignmentOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false);
        final boolean runCoveragePlots = cfgValues.getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true);
        final boolean runCollectBamFileMetrics = cfgValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false);

        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getProject().getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForContext(context);

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = [:]

        for (Sample sample in samples) {

            List<String> availableLibrariesForSample = runtimeService.getLibrariesForSample(sample);
            BamFileGroup mergedBamsPerLibrary = new BamFileGroup();

            // Create per library merged bams
            for (String library in availableLibrariesForSample) {
                BamFileGroup sortedBamFiles = []
                List<LaneFileGroup> rawSequenceGroups = loadLaneFilesForSampleAndLibrary(context, sample, library)
                if (rawSequenceGroups == null || rawSequenceGroups.size() > 0) {
                    for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                        if (runFastQC && !runAlignmentOnly)
                            rawSequenceGroup.calcFastqcForAll();
                        if (runFastQCOnly)
                            continue;

                        BamFile bamFile = rawSequenceGroup.alignAndPairSlim();

                        bamFile.setAsTemporaryFile();  // Bam files created with sai files are only temporary.
                        sortedBamFiles.addFile(bamFile);
                    }
                }

                if (!sortedBamFiles.getFilesInGroup()) continue;

                if (runAlignmentOnly) continue;

                mergedBamsPerLibrary.addFile(sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample));

            }

            // Merge library bams into per sample bams
            BamFile mergedBam = mergedBamsPerLibrary.mergeAndRemoveDuplicatesSlim(sample);
            if (runCollectBamFileMetrics) mergedBam.collectMetrics();

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

    /**
     * Provides a cached method for loading lane files from a sample.
     *
     * @param sample
     * @return
     */
    protected synchronized List<LaneFileGroup> loadLaneFilesForSampleAndLibrary(ExecutionContext context, Sample sample, String library) {
        DataSet dataSet = context.getDataSet();
        if (!foundRawSequenceFileGroups.containsKey(dataSet)) {
            foundRawSequenceFileGroups.put(dataSet, new LinkedHashMap<String, List<LaneFileGroup>>());
        }
        def runtimeService = context.getRuntimeService() as COProjectsRuntimeService
        String sampleID = sample.getName() + "_" + library;
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
                LaneFile copyOfFile = new LaneFile(lf, context);//.getPath(), context, lf.getCreatingJobsResult(), lf.getParentFiles(), lf.getFileStage());
                copyOfFiles.add(copyOfFile);
            }
            copyOfLaneFileGroups.add(new LaneFileGroup(context, lfg.getId(), lfg.getRun(), sample, copyOfFiles));
        }
        return copyOfLaneFileGroups;
    }

    @Override
    public boolean checkExecutability(ExecutionContext context) {
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getProject().getRuntimeService();
        List<Sample> samples = runtimeService.getSamplesForContext(context);
        if (samples.size() == 0)
            return false;

        //Check if at least one file is available. Maybe for two if paired is used...?
        int cnt = 0;
        for (Sample sample : samples) {
            List<LaneFileGroup> laneFileGroups = [];
            for (String lib : sample.getLibraries()) {
                laneFileGroups += loadLaneFilesForSampleAndLibrary(context, sample, lib);
            }

            for (LaneFileGroup lfg : laneFileGroups) {
                cnt += lfg.getFilesInGroup().size();
            }
        }
        return cnt > 0;
    }

}
