package de.dkfz.b080.co.qcworkflow

import de.dkfz.b080.co.common.AlignmentAndQCConfig
import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.b080.co.files.*
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.execution.io.ExecutionService
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import groovy.transform.CompileStatic

/**
 * @author michael
 */
@CompileStatic
class BisulfiteCoreWorkflow extends QCPipeline {

    @Override
    boolean execute(ExecutionContext context) {
        AlignmentAndQCConfig cfg = new AlignmentAndQCConfig(context)
        cfg.setUseOnlyExistingTargetBam(true)

        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.runtimeService
        List<Sample> samples = runtimeService.getSamplesForContext(context)

        BamFileGroup mergedBamFiles = new BamFileGroup()
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = [:]

        for (Sample sample in samples) {
            BamFileGroup mergedBamsPerLibrary = new BamFileGroup()

            // Create per library merged bams
            for (String library in sample.libraries) {
                BamFileGroup sortedBamFiles = []
                List<LaneFileGroup> rawSequenceGroups = runtimeService.
                        loadLaneFilesForSampleAndLibrary(context, sample, library)
                if (rawSequenceGroups == null || rawSequenceGroups.size() > 0) {
                    for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {

                        if (cfg.runFastqc && !cfg.runAlignmentOnly)
                            rawSequenceGroup.calcFastqcForAll()
                        if (cfg.runFastqcOnly)
                            continue


                        BamFile bamFile = rawSequenceGroup.alignAndPairSlim()

                        bamFile.setAsTemporaryFile()  // Bam files created with sai files are only temporary.
                        sortedBamFiles.addFile(bamFile)

                    }
                }

                if (!sortedBamFiles.filesInGroup) continue

                if (cfg.runAlignmentOnly) continue

                BamFile mergedLibraryBam
                if (sample.libraries.size() == 1) {
                    mergedLibraryBam = sortedBamFiles.mergeAndRemoveDuplicatesSlim(sample)
                    if (cfg.runCollectBamFileMetrics) mergedLibraryBam.collectMetrics()

                    if (!coverageTextFilesBySample.containsKey(sample.type))
                        coverageTextFilesBySample.put(sample.type, new CoverageTextFileGroup())
                    coverageTextFilesBySample.get(sample.type).addFile(mergedLibraryBam.calcReadBinsCoverage())

                    mergedBamFiles.addFile(mergedLibraryBam);
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

                mergedBamsPerLibrary.addFile(mergedLibraryBam);
            }

            // Merge library bams into per sample bams
            if(sample.libraries.size() > 1) {
                BamFile mergedBam = mergedBamsPerLibrary.mergeSlim(sample)
                // Unfortunately, due to the way Roddy works, the following call needs to be encapsulated into
                // a method, in order to put library and merged methylation results into different directories.
                // This allows for selection via onMethod="BisulfiteCoreWorkflow.mergedMethylationCallingMeta".
                mergedBam.mergedMethylationCallingMeta()

                if (cfg.runCollectBamFileMetrics) mergedBam.collectMetrics()

                if (!coverageTextFilesBySample.containsKey(sample.type))
                    coverageTextFilesBySample.put(sample.type, new CoverageTextFileGroup())
                coverageTextFilesBySample.get(sample.type).addFile(mergedBam.calcReadBinsCoverage())

                mergedBamFiles.addFile(mergedBam)
            }

        }

        if (!mergedBamFiles.filesInGroup) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.
                    expand("There were no merged bam files available."))
            return false;
        }

        if (cfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() >= 2) {
            coverageTextFilesBySample.get(Sample.SampleType.CONTROL).
                    plotAgainst(coverageTextFilesBySample.get(Sample.SampleType.TUMOR));
        } else if (cfg.runCoveragePlots && coverageTextFilesBySample.keySet().size() == 1) {
            //TODO: Think if this conflicts with plotAgainst on rerun! Maybe missing files are not recognized.
            ((CoverageTextFileGroup) coverageTextFilesBySample.values().toArray()[0]).plot()
        }

        return true
    }


    @Override
    protected boolean checkLaneFiles(ExecutionContext context) {
        boolean returnValue = true
        int cnt = 0;
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        List<Sample> samples = runtimeService.getSamplesForContext(context)
        for (Sample sample : samples) {
            LinkedHashMap<String,LaneFileGroup> laneFileGroups = [:]
            for (String lib : sample.libraries) {
                for (LaneFileGroup group : runtimeService.loadLaneFilesForSampleAndLibrary(context, sample, lib)) {
                    String key = group.run + " " + group.id
                    if (laneFileGroups.containsKey(key)) {
                        context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                                    expand("Duplicate lane identifiers for ${group.run}_${group.id} among libraries " +
                                            "(pid=${context.dataSet}, sample=${sample.name})! Check run and FASTQ names."))
                        returnValue = false
                    } else {
                        laneFileGroups[key] = group
                    }
                    cnt += group.filesInGroup.size()
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
        FileSystemAccessProvider accessProvider = FileSystemAccessProvider.instance
        File cposidx = aqcfg.cytosinePositionIndex
        if (cposidx.toString().equals("")) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("${AlignmentAndQCConfig.CVALUE_CYTOSINE_POSITIONS_INDEX} is not defined!"))
            return false
        } else if (!accessProvider.fileExists(cposidx)
                || !accessProvider.isReadable(cposidx)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("Cytosine position index '${cposidx}' is not accessible!"))
            return false
        } else {
            return true
        }

    }

    protected boolean checkClipIndex(ExecutionContext context) {
        AlignmentAndQCConfig aqcfg = new AlignmentAndQCConfig(context)
        FileSystemAccessProvider accessProvider = FileSystemAccessProvider.instance
        if (!aqcfg.clipIndex.toString().startsWith('$' + ExecutionService.RODDY_CVALUE_DIRECTORY_EXECUTION)
                && !accessProvider.isReadable(aqcfg.clipIndex)) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_SETUP_INVALID.
                    expand("Clip index '${aqcfg.clipIndex}' is not accessible!"))
            return false
        } else {
            return true
        }
    }

    @Override
    boolean checkExecutability(ExecutionContext context) {
        boolean result = super.checkExecutability(context)
        result &= checkCytosinePositionIndex(context)
        result &= checkClipIndex(context)
        return result
    }

}
