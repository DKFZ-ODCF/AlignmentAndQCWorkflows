package de.dkfz.b080.co.common;

import de.dkfz.b080.co.files.*
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.core.DataSet
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.core.ExecutionContextLevel
import de.dkfz.roddy.core.ExecutionContextSubLevel
import de.dkfz.roddy.core.ProcessingFlag
import de.dkfz.roddy.execution.io.ExecutionResult;
import de.dkfz.roddy.execution.io.ExecutionService;
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.tools.LoggerWrapper
import groovy.transform.CompileStatic

import java.util.function.Consumer;


@CompileStatic
class COProjectsRuntimeService extends BasicCOProjectsRuntimeService {
    private static LoggerWrapper logger = LoggerWrapper.getLogger(BasicCOProjectsRuntimeService.class.getName());

    protected static void getFileCompression(ExecutionContext run, List<LaneFile> allLaneFiles) {
        ExecutionService es = ExecutionService.getInstance();
        StringBuilder command = new StringBuilder();
        for (LaneFile lf : allLaneFiles) {
            if (lf.getDecompressionString() == null) {
                command << "&& file -bL ${lf.path.absolutePath}";
//            } else if (lf.getDecompressionString().equals("setgid")) {
//                command << "&& file -bL ${lf.path.absolutePath} | cut -d' ' -f2 ";
            } else {
                command << "&& echo 0 ";
            }
        }

        //TODO Throw or post error that no lane files are available

        ExecutionResult er = es.execute(command[2..-1].toString());
        for (int i = 0; i < er.resultLines.size(); i++) {
            //Look at Matthias ExomePipeline extension. The code is mainly taken from there.
            String dCString = null;
            String dRString = "gzip -c"; //POSSIBLE-ERROR: Not zipper set for 'cat'/ASCII
            String type = er.resultLines[i];
            if (type.startsWith("setgid "))
                type = type[7..-1].trim().split(" ")[0];
            if (type == "0") continue;
            if (type == "gzip") {
                dCString = "gunzip -c";
                dRString = "gzip -c"
            } else if (type == "bzip2") {
                dCString = "bunzip2 -c -k";
                dRString = "bzip2 -c -k";
            } else if (type == "ASCII" || type == "ASCII text") {
                dCString = "cat";
                dRString = "head -n E";
            }
            allLaneFiles[i].setDecompressionString(dCString);
            allLaneFiles[i].setRecompressionString(dRString);
        }
    }

    public static void determineFileDecoder(ExecutionContext context, List<LaneFileGroup> laneFileGroupList) {
        List<LaneFile> allLaneFiles = [];
        for (LaneFileGroup lfg : laneFileGroupList) {
            List<LaneFile> lflist = lfg.getFilesInGroup();
            for (LaneFile lf : lflist) {
                allLaneFiles << lf;
            }
        }
        getFileCompression(context, allLaneFiles);
    }

    public List<LaneFileGroup> getLanesForSampleAndLibrary(ExecutionContext context, Sample sample, String library) {
        return getLanesForSample(context, sample, library);
    }

    public List<LaneFileGroup> getLanesForSample(ExecutionContext context, Sample sample, String libraryID = null) {
        AlignmentAndQCConfig config = new AlignmentAndQCConfig(context);

        ProcessingFlag flag = context.setProcessingFlag(ProcessingFlag.STORE_FILES);

        List<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>()
        if (config.getMetadataTableIsSet()) {
            laneFiles = getLaneFileGroupsFromInputTable(context, sample, libraryID)
        } else if (config.getFastqFileListIsSet()) {
            laneFiles = getLaneFileGroupsFromFastqList(context, sample, libraryID)
        } else {
            laneFiles = getLaneFileGroupsFromFilesystem(context, sample, libraryID)
        }
        if (laneFiles.size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no lane files available for sample ${sample.getName()}"));
        } else {
            determineFileDecoder(context, laneFiles);
        }
        context.setProcessingFlag(flag);
        return laneFiles;
    }

    /**
     * This entry is used for caching purposes.
     */
    protected Map<DataSet, Map<String, List<LaneFileGroup>>> foundRawSequenceFileGroups = new LinkedHashMap<>();

    public synchronized List<LaneFileGroup> loadLaneFilesForSample(ExecutionContext context, Sample sample) {
        DataSet dataSet = context.getDataSet();
        if (!foundRawSequenceFileGroups.containsKey(dataSet)) {
            foundRawSequenceFileGroups.put(dataSet, new LinkedHashMap<String, List<LaneFileGroup>>());
        }
        COProjectsRuntimeService runtimeService = (COProjectsRuntimeService) context.getRuntimeService();
        String sampleID = sample.getName();
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
                LaneFile copyOfFile = new LaneFile(lf, context);//lf.getPath(), context, lf.getCreatingJobsResult(), lf.getParentFiles(), lf.getFileStage());
                copyOfFiles.add(copyOfFile);
            }
            copyOfLaneFileGroups.add(new LaneFileGroup(context, lfg.getId(), lfg.getRun(), sample, copyOfFiles));
        }
        return copyOfLaneFileGroups;
    }


    /**
     * Provides a cached method for loading lane files from a sample.
     *
     * @param sample
     * @return
     */
    public synchronized List<LaneFileGroup> loadLaneFilesForSampleAndLibrary(ExecutionContext context, Sample sample, String library) {
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

    public List<LaneFileGroup> getLaneFileGroupsFromInputTable(ExecutionContext context, Sample sample, String libraryID = null) {
        MetadataTable inputTable = getMetadataTable(context).subsetBySample(sample.name)
        if (libraryID) inputTable = inputTable.subsetByLibrary(libraryID)
        List<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>()
        for(String runID : inputTable.listRunIDs()) {
            List<File> fastqFilesForRun = inputTable.subsetByRun(runID).listFiles()
            List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, runID, fastqFilesForRun);
            laneFiles.addAll(bundleFiles)
        }
        return laneFiles
    }

    protected static int indexOfPathElement(String pathnamePattern, String element) {
        int index = pathnamePattern.split(StringConstants.SPLIT_SLASH).findIndexOf { it -> it == element }
        if (index < 0) {
            throw new RuntimeException("Couldn't match '${element}' in '${pathnamePattern}")
        }
        return index
    }

    public List<LaneFileGroup> getLaneFileGroupsFromFastqList(ExecutionContext context, Sample sample, String libraryID) {
        COConfig coConfig = new COConfig(context)
        List<File> fastqFiles = coConfig.getFastqList().collect { String it -> new File(it); }
        def sequenceDirectory = coConfig.getSequenceDirectory()
        int indexOfSampleID = indexOfPathElement(sequenceDirectory, '${sample}')
        int indexOfRunID = indexOfPathElement(sequenceDirectory, '${run}')
        List<File> listOfFastqFilesForSample = fastqFiles.findAll {
            it.absolutePath.split(StringConstants.SPLIT_SLASH)[indexOfSampleID] == sample.getName()
        } as List<File>
        List<String> listOfRunIDs = listOfFastqFilesForSample.collect {
            it.absolutePath.split(StringConstants.SPLIT_SLASH)[indexOfRunID]
        }.unique() as List<String>
        List<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>()
        for(String runID : listOfRunIDs) {
            List<File> fastqFilesForRun = listOfFastqFilesForSample.findAll {
                it.absolutePath.split(StringConstants.SPLIT_SLASH)[indexOfRunID] == runID
            } as List<File>
            List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, runID, fastqFilesForRun);
            laneFiles.addAll(bundleFiles)
        }
        return laneFiles
    }

    public List<LaneFileGroup> getLaneFileGroupsFromFilesystem(ExecutionContext context, Sample sample, String libraryID) {
        File sampleDirectory = getSampleDirectory(context, sample, libraryID);

        logger.postAlwaysInfo("Searching for lane files in directory ${sampleDirectory}")
        List<File> runsForSample = FileSystemAccessProvider.getInstance().listDirectoriesInDirectory(sampleDirectory);
        LinkedList<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>()
        for(File run : runsForSample) {
            File sequenceDirectory = getSequenceDirectory(context, sample, run.getName(), libraryID);
            logger.postSometimesInfo("\tChecking for run ${run.name} in sequence dir: ${sequenceDirectory}")
            if (!FileSystemAccessProvider.getInstance().checkDirectory(sequenceDirectory, context, false)) // Skip directories which do not exist
                continue;
            List<File> files = FileSystemAccessProvider.getInstance().listFilesInDirectory(sequenceDirectory);
            if (files.size() == 0)
                logger.postAlwaysInfo("\t There were no lane files in directory ${sequenceDirectory}")
            //Find file bundles
            List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, run.getName(), files);
            logger.postSometimesInfo("\tFound ${bundleFiles.size()} lane file groups in sequence directory.")
            laneFiles.addAll(bundleFiles)
        }
        return laneFiles
    }

    public BamFileGroup getPairedBamFilesForDataSet(ExecutionContext context, Sample sample) {
        File alignmentDirectory = getAlignmentDirectory(context);
        final String pairedBamSuffix = context.getConfiguration().getConfigurationValues().get("pairedBamSuffix", "paired.bam.sorted.bam")
        //TODO Create constants
        List<String> filters = ["${sample.getName()}*${pairedBamSuffix}".toString()]
        List<File> pairedBamPaths = FileSystemAccessProvider.getInstance().listFilesInDirectory(alignmentDirectory, filters);

        int laneID = 0;
        List<BamFile> bamFiles = pairedBamPaths.collect({
            File f ->
                laneID++;
                String name = f.getName();
                String[] split = name.split(StringConstants.SPLIT_UNDERSCORE);
                String sampleName = split[0];
                int runIndex = 1;
                if (split[1].isInteger()) {
                    runIndex = 2;
                    sampleName = split[0..1].join(StringConstants.UNDERSCORE);
                }
                RunID run = new RunID(split[runIndex..-2].join(StringConstants.UNDERSCORE))
                LaneID lane = new LaneID(String.format("L%03d", laneID))


                BamFile bamFile =
                        COBaseFile.constructSourceFile(
                                BamFile,
                                f,
                                context,
                                new COFileStageSettings(
                                        lane,
                                        run,
                                        (LibraryID) null,
                                        sample,
                                        context.dataSet)
                        ) as BamFile
                return bamFile;
        })
        BamFileGroup bamFileGroup = new BamFileGroup(bamFiles);
        logger.info("Found ${bamFileGroup.getFilesInGroup().size()} paired bam files for sample ${sample.getName()}");
        return bamFileGroup;
    }

    public BamFileGroup getMergedBamFilesForDataSet(ExecutionContext context) {

    }

    public BamFileGroup getMergedBamFilesForDataSet(ExecutionContext context, DataSet dataSet) {

        ProcessingFlag flag = context.setProcessingFlag(ProcessingFlag.STORE_NOTHING);

        BamFileGroup mergedBamFiles = new BamFileGroup();
        Map<Sample.SampleType, CoverageTextFileGroup> coverageTextFilesBySample = new LinkedHashMap<>();

        //Set level to test, set back later.
        ExecutionContextLevel executionContextLevel = context.getExecutionContextLevel();
        ExecutionContextSubLevel detailedExecutionContextLevel = context.getDetailedExecutionContextLevel();
        context.setExecutionContextLevel(ExecutionContextLevel.QUERY_STATUS);
//        context.setExecutionContextLevel(ExecutionContextLevel.QUERY_STATUS);
        getSamplesForContext(context).parallelStream().forEach(new Consumer<Sample>() {
            @Override
            void accept(Sample sample) {
                List<LaneFileGroup> rawSequenceGroups = getLanesForSample(context, sample);
                BamFileGroup sortedBamFiles = new BamFileGroup();

                for (LaneFileGroup rawSequenceGroup : rawSequenceGroups) {
                    rawSequenceGroup.alignAll();
                    BamFile bamFile = rawSequenceGroup.getAllAlignedFiles().pairAndSort();
                    sortedBamFiles.addFile(bamFile);
                }
                context.setProcessingFlag(ProcessingFlag.STORE_FILES);
                mergedBamFiles.addFile(sortedBamFiles.mergeAndRemoveDuplicates());
                context.setProcessingFlag(ProcessingFlag.STORE_NOTHING);
            }
        });
        context.setExecutionContextLevel(executionContextLevel);
        context.setProcessingFlag(flag);
        return mergedBamFiles;
    }

}
