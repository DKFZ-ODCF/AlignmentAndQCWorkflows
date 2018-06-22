/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.common

import de.dkfz.b080.co.files.*
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.core.DataSet
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.core.ProcessingFlag
import de.dkfz.roddy.execution.io.ExecutionResult
import de.dkfz.roddy.execution.io.ExecutionService
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.tools.LoggerWrapper
import de.dkfz.roddy.tools.RoddyIOHelperMethods

/**
 * Created by heinold on 15.01.16.
 */
@groovy.transform.CompileStatic
class COProjectsRuntimeService extends BasicCOProjectsRuntimeService {
    private static LoggerWrapper logger = LoggerWrapper.getLogger(BasicCOProjectsRuntimeService.class.getName())

    protected static void setFileCompressionInLaneFiles(ExecutionContext run, List<LaneFile> allLaneFiles) {
        ExecutionService es = ExecutionService.getInstance()
        StringBuilder command = new StringBuilder()
        for (LaneFile lf : allLaneFiles) {
            if (lf.getDecompressionString() == null) {
                command << "&& file -bL ${lf.path.absolutePath}"
//            } else if (lf.getDecompressionString().equals("setgid")) {
//                command << "&& file -bL ${lf.path.absolutePath} | cut -d' ' -f2 "
            } else {
                command << "&& echo 0 "
            }
        }

        //TODO Throw or post error that no lane files are available

        ExecutionResult er = es.execute(command[2..-1].toString())
        for (int i = 0; i < er.resultLines.size(); i++) {
            //Look at Matthias ExomePipeline extension. The code is mainly taken from there.
            String dCString = null
            String dRString = "gzip -c" //POSSIBLE-ERROR: Not zipper set for 'cat'/ASCII
            String type = er.resultLines[i]
            if (type.startsWith("setgid "))
                type = type[7..-1].trim().split(" ")[0]
            if (type == "0") continue
            if (type == "gzip") {
                dCString = "gunzip -c"
                dRString = "gzip -c"
            } else if (type == "bzip2") {
                dCString = "bunzip2 -c -k"
                dRString = "bzip2 -c -k"
            } else if (type == "ASCII" || type == "ASCII text") {
                dCString = "cat"
                dRString = "head -n E"
            }
            allLaneFiles[i].setDecompressionString(dCString)
            allLaneFiles[i].setRecompressionString(dRString)
        }
    }

    static void setFileCompressionInLaneFileGroups(ExecutionContext context, List<LaneFileGroup> laneFileGroupList) {
        List<LaneFile> allLaneFiles = []
        for (LaneFileGroup lfg : laneFileGroupList) {
            List<LaneFile> lflist = lfg.getFilesInGroup()
            for (LaneFile lf : lflist) {
                allLaneFiles << lf
            }
        }
        setFileCompressionInLaneFiles(context, allLaneFiles)
    }

    List<LaneFileGroup> getLanesForSample(ExecutionContext context, Sample sample, String libraryID = null) {
        AlignmentAndQCConfig config = new AlignmentAndQCConfig(context)

        ProcessingFlag flag = context.setProcessingFlag(ProcessingFlag.STORE_FILES)

        List<LaneFileGroup> laneFiles
        if (config.getMetadataTableIsSet()) {
            laneFiles = getLaneFileGroupsFromInputTable(context, sample, libraryID)
        } else if (config.getFastqFileListIsSet()) {
            laneFiles = getLaneFileGroupsFromFastqList(context, sample, libraryID)
        } else {
            laneFiles = getLaneFileGroupsFromFilesystem(context, sample, libraryID)
        }
        if (laneFiles.size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no lane files available for sample ${sample.getName()}"))
        } else {
            setFileCompressionInLaneFileGroups(context, laneFiles)
        }
        context.setProcessingFlag(flag)
        return laneFiles
    }

    /**
     * This entry is used for caching purposes.
     */
    protected Map<DataSet, Map<String, List<LaneFileGroup>>> foundRawSequenceFileGroups = new LinkedHashMap<>()

    synchronized List<LaneFileGroup> loadLaneFilesForSample(ExecutionContext context, Sample sample) {
        DataSet dataSet = context.getDataSet()
        if (!foundRawSequenceFileGroups.containsKey(dataSet)) {
            foundRawSequenceFileGroups.put(dataSet, new LinkedHashMap<String, List<LaneFileGroup>>())
        }
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService()
        String sampleID = sample.getName()
        Map<String, List<LaneFileGroup>> mapForDataSet = foundRawSequenceFileGroups.get(dataSet)
        if (!mapForDataSet.containsKey(sampleID)) {
            List<LaneFileGroup> laneFileGroups = runtimeService.getLanesForSample(context, sample)
            mapForDataSet.put(sampleID, laneFileGroups)
        }

        List<LaneFileGroup> laneFileGroups = mapForDataSet.get(sampleID)
        List<LaneFileGroup> copyOfLaneFileGroups = new LinkedList<LaneFileGroup>()
        for (LaneFileGroup lfg : laneFileGroups) {
            List<LaneFile> copyOfFiles = new LinkedList<>()
            for (LaneFile lf : lfg.getFilesInGroup()) {
                LaneFile copyOfFile = new LaneFile(lf, context)//lf.getPath(), context, lf.getCreatingJobsResult(), lf.getParentFiles(), lf.getFileStage())
                copyOfFiles.add(copyOfFile)
            }
            copyOfLaneFileGroups.add(new LaneFileGroup(context, lfg.getId(), lfg.getRun(), sample, copyOfFiles))
        }
        return copyOfLaneFileGroups
    }


    /**
     * Provides a cached method for loading lane files from a sample.
     *
     * @param sample
     * @return
     */
    synchronized List<LaneFileGroup> loadLaneFilesForSampleAndLibrary(ExecutionContext context, Sample sample, String library) {
        DataSet dataSet = context.getDataSet()
        if (!foundRawSequenceFileGroups.containsKey(dataSet)) {
            foundRawSequenceFileGroups.put(dataSet, new LinkedHashMap<String, List<LaneFileGroup>>())
        }
        def runtimeService = context.getRuntimeService() as COProjectsRuntimeService
        String sampleID = sample.getName() + "_" + library
        Map<String, List<LaneFileGroup>> mapForDataSet = foundRawSequenceFileGroups.get(dataSet)
        if (!mapForDataSet.containsKey(sampleID)) {
            List<LaneFileGroup> laneFileGroups = runtimeService.getLanesForSample(context, sample, library)
            mapForDataSet.put(sampleID, laneFileGroups)
        }

        List<LaneFileGroup> laneFileGroups = mapForDataSet.get(sampleID)
        List<LaneFileGroup> copyOfLaneFileGroups = new LinkedList<LaneFileGroup>()
        for (LaneFileGroup lfg : laneFileGroups) {
            List<LaneFile> copyOfFiles = new LinkedList<>()
            for (LaneFile lf : lfg.getFilesInGroup()) {
                LaneFile copyOfFile = new LaneFile(lf, context)
                copyOfFiles.add(copyOfFile)
            }
            copyOfLaneFileGroups.add(new LaneFileGroup(context, lfg.getId(), lfg.getRun(), sample, copyOfFiles))
        }
        return copyOfLaneFileGroups
    }

    List<String> listSampleNames(MetadataTable inputTable) {
        return inputTable.listColumn(COConstants.INPUT_TABLE_SAMPLECOL_NAME).unique()
    }

    List<String> listRunIds(MetadataTable inputTable) {
        return inputTable.listColumn(COConstants.INPUT_TABLE_RUNCOL_NAME).unique()
    }

    List<String> listLibraries(MetadataTable inputTable) {
        return inputTable.listColumn(COConstants.INPUT_TABLE_MARKCOL_NAME).unique()
    }


    List<LaneFileGroup> getLaneFileGroupsFromInputTable(ExecutionContext context, Sample sample, String libraryID = null) {
        MetadataTable inputTable = metadataAccessor.getMetadataTable(context).subsetBySample(sample.name)
        if (libraryID) inputTable = inputTable.subsetByLibrary(libraryID)
        List<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>()
        for(String runID : listRunIds(inputTable)) {
            List<File> fastqFilesForRun = inputTable.subsetByRun(runID).listFiles()
            List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, runID, fastqFilesForRun)
            laneFiles.addAll(bundleFiles)
        }
        return laneFiles
    }

    List<LaneFileGroup> getLaneFileGroupsFromFastqList(ExecutionContext context, Sample sample, String libraryID) {
        COConfig coConfig = new COConfig(context)
        List<File> fastqFiles = coConfig.getFastqList().collect { String it -> new File(it) }
        def sequenceDirectory = coConfig.getSequenceDirectory()
        Integer indexOfSampleID = RoddyIOHelperMethods.findComponentIndexInPath(sequenceDirectory, '${sample}').
                orElseThrow { new RuntimeException("Could not find '\${sample}' in ${COConstants.CVALUE_SEQUENCE_DIRECTORY}='${sequenceDirectory}'") }
        Integer indexOfRunID = RoddyIOHelperMethods.findComponentIndexInPath(sequenceDirectory, '${run}').
                orElseThrow { new RuntimeException("Could not finde '\${run}' in ${COConstants.CVALUE_SEQUENCE_DIRECTORY}='${sequenceDirectory}'") }
        List<File> listOfFastqFilesForSample = fastqFiles.findAll {
            RoddyIOHelperMethods.splitPathname(it.absolutePath)[indexOfSampleID] == sample.getName()
        } as List<File>
        List<String> listOfRunIDs = listOfFastqFilesForSample.collect {
            RoddyIOHelperMethods.splitPathname(it.absolutePath)[indexOfRunID]
        }.unique() as List<String>
        List<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>()
        for(String runID : listOfRunIDs) {
            List<File> fastqFilesForRun = listOfFastqFilesForSample.findAll {
                RoddyIOHelperMethods.splitPathname(it.absolutePath)[indexOfRunID] == runID
            } as List<File>
            List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, runID, fastqFilesForRun)
            laneFiles.addAll(bundleFiles)
        }
        return laneFiles
    }

    List<LaneFileGroup> getLaneFileGroupsFromFilesystem(ExecutionContext context, Sample sample, String libraryID) {
        File sampleDirectory = getSampleDirectory(context, sample, libraryID)

        logger.postAlwaysInfo("Searching for lane files in directory ${sampleDirectory}")
        List<File> runsForSample = FileSystemAccessProvider.getInstance().listDirectoriesInDirectory(sampleDirectory)
        LinkedList<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>()
        for(File run : runsForSample) {
            File sequenceDirectory = getSequenceDirectory(context, sample, run.getName(), libraryID)
            logger.postSometimesInfo("\tChecking for run ${run.name} in sequence dir: ${sequenceDirectory}")
            if (!FileSystemAccessProvider.getInstance().checkDirectory(sequenceDirectory, context, false)) // Skip directories which do not exist
                continue
            List<File> files = FileSystemAccessProvider.getInstance().listFilesInDirectory(sequenceDirectory)
            if (files.size() == 0)
                logger.postAlwaysInfo("\t There were no lane files in directory ${sequenceDirectory}")
            //Find file bundles
            List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, run.getName(), files)
            logger.postSometimesInfo("\tFound ${bundleFiles.size()} lane file groups in sequence directory.")
            laneFiles.addAll(bundleFiles)
        }
        return laneFiles
    }

    BamFileGroup getPairedBamFilesForDataSet(ExecutionContext context, Sample sample) {
        File alignmentDirectory = getAlignmentDirectory(context)
        final String pairedBamSuffix = context.getConfiguration().getConfigurationValues().get("pairedBamSuffix", "paired.bam.sorted.bam")
        //TODO Create constants
        List<String> filters = ["${sample.getName()}*${pairedBamSuffix}".toString()]
        List<File> pairedBamPaths = FileSystemAccessProvider.getInstance().listFilesInDirectory(alignmentDirectory, filters)

        int laneID = 0
        List<BamFile> bamFiles = pairedBamPaths.collect({
            File f ->
                laneID++
                String name = f.getName()
                String[] split = name.split(StringConstants.SPLIT_UNDERSCORE)
                String sampleName = split[0]
                int runIndex = 1
                if (split[1].isInteger()) {
                    runIndex = 2
                    sampleName = split[0..1].join(StringConstants.UNDERSCORE)
                }
                String run = split[runIndex..-2].join(StringConstants.UNDERSCORE)
                String lane = String.format("L%03d", laneID)

                BamFile bamFile = COBaseFile.constructSourceFile(BamFile, f, context,
                        new COFileStageSettings(new LaneID(lane), new RunID(run), null, sample, context.getDataSet())) as BamFile
                return bamFile
        })
        logger.info("Found ${bamFiles.size()} paired bam files for sample ${sample.getName()}")
        return new BamFileGroup(bamFiles)
    }
    
}
