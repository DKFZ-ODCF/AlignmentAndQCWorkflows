package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.common.AlignmentAndQCConfig
import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.BamFile
import de.dkfz.b080.co.files.BamFileGroup
import de.dkfz.b080.co.files.CoverageTextFileGroup
import de.dkfz.b080.co.files.LaneFileGroup
import de.dkfz.b080.co.files.Sample
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider

/**
 * @author michael
 */
@groovy.transform.CompileStatic
class BisulfiteCoreWorkflow extends QCPipeline {

    @Override
    boolean execute(ExecutionContext context) {
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        aqcfg.extractSamplesFromOutputFiles = false

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)

        BamFileGroup mergedBamFiles = new BamFileGroup()
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = [:]

        for (Sample sample in samples) {

            List<String> availableLibrariesForSample = sample.getLibraries()
            BamFileGroup mergedBamsPerLibrary = new BamFileGroup()

            // Create per library merged bams
            for (String library in availableLibrariesForSample) {
                BamFileGroup sortedBamFiles = []
                List<LaneFileGroup> rawSequenceGroups = runtimeService.loadLaneFilesForSampleAndLibrary(context, sample, library)
                if (rawSequenceGroups == null || rawSequenceGroups.size() > 0) {
                    for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                        LinkedHashMap<String, String> libraryParameters = ["library": library, "LIBRARY": library]

                        if (aqcfg.runFastQC && !aqcfg.runAlignmentOnly)
                            rawSequenceGroup.calcFastqcForAll(libraryParameters)
                        if (aqcfg.runFastQCOnly)
                            continue

                        BamFile bamFile = rawSequenceGroup.alignAndPairSlim(libraryParameters)

                        bamFile.setAsTemporaryFile()  // Bam files created with sai files are only temporary.
                        sortedBamFiles.addFile(bamFile)

                    }
                }

                if (!sortedBamFiles.getFilesInGroup()) continue

                if (aqcfg.runAlignmentOnly) continue

                BamFile mergedLibraryBam
                if (availableLibrariesForSample.size() == 1) {
                    mergedLibraryBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample)
                    if (aqcfg.runCollectBamFileMetrics) mergedLibraryBam.collectMetrics()

                    Sample.SampleType sampleType = sample.getType()
                    if (!coverageTextFilesBySample.containsKey(sampleType))
                        coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup())
                    coverageTextFilesBySample.get(sampleType).addFile(mergedLibraryBam.calcReadBinsCoverage())

                    mergedBamFiles.addFile(mergedLibraryBam)
                    // Unfortunately, due to the way Roddy works, the following call needs to be encapsulated into
                    // a method, in order to put library and merged methylation results into different directories.
                    // This allows for selection via onMethod="BisulfiteCoreWorkflow.mergedMethylationCallingMeta".
                    mergedLibraryBam.mergedMethylationCallingMeta()
                }
                else {
                    mergedLibraryBam = sortedBamFiles.mergeAndRemoveDuplicatesSlimWithLibrary(sample, library)
                    // Unfortunately, due to the way Roddy works, the following call needs to be encapsulated into
                    // a method, in order to put library and merged methylation results into different directories.
                    // This allows for selection via onMethod="BisulfiteCoreWorkflow.libraryMethylationCallingMeta".
                    mergedLibraryBam.libraryMethylationCallingMeta()
                }

                mergedBamsPerLibrary.addFile(mergedLibraryBam)
            }

            // Merge library bams into per sample bams
            if(availableLibrariesForSample.size() > 1) {
                BamFile mergedBam = mergedBamsPerLibrary.mergeSlim(sample)
                // Unfortunately, due to the way Roddy works, the following call needs to be encapsulated into
                // a method, in order to put library and merged methylation results into different directories.
                // This allows for selection via onMethod="BisulfiteCoreWorkflow.mergedMethylationCallingMeta".
                mergedBam.mergedMethylationCallingMeta()

                if (aqcfg.runCollectBamFileMetrics) mergedBam.collectMetrics()

                Sample.SampleType sampleType = sample.getType()
                if (!coverageTextFilesBySample.containsKey(sampleType))
                    coverageTextFilesBySample.put(sampleType, new CoverageTextFileGroup())
                coverageTextFilesBySample.get(sampleType).addFile(mergedBam.calcReadBinsCoverage())

                mergedBamFiles.addFile(mergedBam)
            }

        }

        if (!mergedBamFiles.getFilesInGroup()) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no merged bam files available."))
            return false
        }

        if (aqcfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR))
        } else if (aqcfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot()
        }

        return true
    }


    @Override
    protected boolean checkLaneFiles(ExecutionContext context) {
        boolean returnValue = true
        int cnt = 0
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        for (Sample sample : samples) {
            LinkedHashMap<String,LaneFileGroup> laneFileGroups = [:]
            for (String lib : sample.getLibraries()) {
                for (LaneFileGroup group : runtimeService.loadLaneFilesForSampleAndLibrary(context, sample, lib)) {
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

    protected boolean checkCytosinePositionIndex(ExecutionContext context) {
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        FileSystemAccessProvider accessProvider = FileSystemAccessProvider.getInstance()
        File cposidx = aqcfg.getCytosinePositionIndex()
        if (cposidx.toString().equals("")) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("${AlignmentAndQCConfig.CVALUE_CYTOSINE_POSITIONS_INDEX} is not defined!"))
            return false
        } else if (!accessProvider.fileExists(cposidx)
                || !accessProvider.isReadable(cposidx)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.expand("Cytosine position index '${cposidx}' is not accessible!"))
            return false
        } else {
            return true
        }

    }

    @Override
    boolean checkExecutability(ExecutionContext context) {
        boolean result = super.checkExecutability(context)
        result &= checkCytosinePositionIndex(context)
        return result
    }

}
