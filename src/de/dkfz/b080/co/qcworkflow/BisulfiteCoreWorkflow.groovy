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
        //cfgValues.put(COConstants.FLAG_USE_ACCELERATED_HARDWARE, "false", "boolean"); //Disable accelerated hardware usage for testing

        // Run flags
        final boolean runFastQCOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC_ONLY, false)
        final boolean runFastQC = cfgValues.getBoolean(COConstants.FLAG_RUN_FASTQC, true)
        final boolean runAlignmentOnly = cfgValues.getBoolean(COConstants.FLAG_RUN_ALIGNMENT_ONLY, false);
        final boolean runCoveragePlots = cfgValues.getBoolean(COConstants.FLAG_RUN_COVERAGE_PLOTS, true);
        final boolean runCollectBamFileMetrics = cfgValues.getBoolean(COConstants.FLAG_RUN_COLLECT_BAMFILE_METRICS, false);

        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService();

        List<Sample> samples = runtimeService.getSamplesForContext(context);

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = [:]

        for (Sample sample in samples) {
            List<String> availableLibrariesForSample = sample.getLibraries();
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
                        //rawSequenceGroup.alignAndPairSlim()


                        bamFile.setAsTemporaryFile();  // Bam files created with sai files are only temporary.
                        sortedBamFiles.addFile(bamFile);

                    }
                }

                if (!sortedBamFiles.getFilesInGroup()) continue;

                if (runAlignmentOnly) continue;

                BamFile mergedLibraryBam;
                if (availableLibrariesForSample.size() == 1) {
                    mergedLibraryBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample);
                    if (runCollectBamFileMetrics) mergedLibraryBam.collectMetrics();

                    Sample.SampleType sampleType = sample.getType();
                    if (!coverageTextFilesBySample.containsKey(sampleType))
                        coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup());
                    coverageTextFilesBySample.get(sampleType).addFile(mergedLibraryBam.calcReadBinsCoverage());

                    mergedBamFiles.addFile(mergedLibraryBam);
                }
                else {
                    mergedLibraryBam = sortedBamFiles.mergeAndRemoveDuplicatesSlimWithLibrary(sample, library)
                }

                mergedBamsPerLibrary.addFile(mergedLibraryBam);
                // Unfortunately, due to the way Roddy works, the following call needs to be encapsulated into
                // a method, in order to put library and merged methylation results into different directories.
                // This allows for selection via onMethod="BisulfiteCoreWorkflow.libraryMethylationCallingMeta".
                mergedLibraryBam.libraryMethylationCallingMeta()
            }

            // Merge library bams into per sample bams
            if(availableLibrariesForSample.size() > 1) {
                BamFile mergedBam = mergedBamsPerLibrary.mergeSlim(sample);
                // Unfortunately, due to the way Roddy works, the following call needs to be encapsulated into
                // a method, in order to put library and merged methylation results into different directories.
                // This allows for selection via onMethod="BisulfiteCoreWorkflow.mergedMethylationCallingMeta".
                mergedBam.mergedMethylationCallingMeta()

                if (runCollectBamFileMetrics) mergedBam.collectMetrics();

                Sample.SampleType sampleType = sample.getType();
                if (!coverageTextFilesBySample.containsKey(sampleType))
                    coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup());
                coverageTextFilesBySample.get(sampleType).addFile(mergedBam.calcReadBinsCoverage());

                mergedBamFiles.addFile(mergedBam);
            }

        }

        if (!mergedBamFiles.getFilesInGroup()) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."));
            return false;
        }

        if (runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (runCoveragePlots && coverageTextFilesBySample.keySet().size() == 1) {
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
            List<LaneFileGroup> laneFileGroups = runtimeService.getLanesForSample(context, sample, library);
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
    protected boolean checkLaneFiles(ExecutionContext context) {
        boolean returnValue = true
        int cnt = 0;
        BasicCOProjectsRuntimeService runtimeService = (BasicCOProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        for (Sample sample : samples) {
            LinkedHashMap<String,LaneFileGroup> laneFileGroups = [:]
            for (String lib : sample.getLibraries()) {
                for (LaneFileGroup group : loadLaneFilesForSampleAndLibrary(context, sample, lib)) {
                    String key = group.getRun() + " " + group.getId()
                    if (laneFileGroups.containsKey(key)) {
                        context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                                    expand("Duplicate lane identifiers for ${group.getRun()}_${group.getId()} among libraries " +
                                            "(pid=${context.getDataSet()}, sample=${sample.getName()})! Check run and FASTQ names."))
                        returnValue = false
                    } else {
                        laneFileGroups[key] = group
                    }
                    cnt += group.getFilesInGroup().size()
                }
            }
        }
        if (cnt <= 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                        expand("No lane files found for PID ${context.getDataSet()}!"))
            returnValue = false
        }
        return returnValue
    }

}
