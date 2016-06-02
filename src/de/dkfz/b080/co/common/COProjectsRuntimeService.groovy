package de.dkfz.b080.co.common;

import de.dkfz.b080.co.files.*
import de.dkfz.roddy.Roddy;
import de.dkfz.roddy.StringConstants
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.core.ExecutionContextLevel
import de.dkfz.roddy.core.ExecutionContextSubLevel
import de.dkfz.roddy.core.ProcessingFlag
import de.dkfz.roddy.execution.io.ExecutionResult;
import de.dkfz.roddy.execution.io.ExecutionService;
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.tools.LoggerWrapper

import java.util.function.Consumer;

/**
 * Created by heinold on 15.01.16.
 */
@groovy.transform.CompileStatic
public class COProjectsRuntimeService extends BasicCOProjectsRuntimeService {
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

    public List getLanesForSample(ExecutionContext context, Sample sample, String libraryID = null) {
        COConfig coConfig = new COConfig(context);

        ProcessingFlag flag = context.setProcessingFlag(ProcessingFlag.STORE_FILES);
        List<LaneFileGroup> laneFiles = new LinkedList<LaneFileGroup>();

        def configurationValues = context.getConfiguration().getConfigurationValues()
        boolean getLanesFromFastqList = configurationValues.getString("fastq_list", "");
        boolean getLanesFromInputTable = Roddy.isMetadataCLOptionSet()

        if (getLanesFromInputTable) {
            MetadataTable inputTable = getMetadataTable(context).subsetBySample(sample.name)
            if (libraryID) inputTable = inputTable.subsetByLibrary(libraryID)

            inputTable.listRunIDs().each {
                String runID ->
                    List<File> fastqFilesForRun = inputTable.subsetByColumn(COConstants.INPUT_TABLE_RUNCOL_NAME, runID).listFiles()
                    List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, runID, fastqFilesForRun);
                    laneFiles += bundleFiles
            }

        } else if (getLanesFromFastqList) {
            // If fastq_list was set via command line or via config file.
            List<File> fastqFiles = configurationValues.getString("fastq_list").split(StringConstants.SPLIT_SEMICOLON).collect { String it -> new File(it); };
            def sequenceDirectory = configurationValues.get(COConstants.CVALUE_SEQUENCE_DIRECTORY).toFile(context).getAbsolutePath();
            int indexOfSampleID = sequenceDirectory.split(StringConstants.SPLIT_SLASH).findIndexOf { it -> it == '${sample}' }
            int indexOfRunID = sequenceDirectory.split(StringConstants.SPLIT_SLASH).findIndexOf { it -> it == '${run}' }

            List<File> listOfFastqFilesForSample = fastqFiles.findAll { it.absolutePath.split(StringConstants.SPLIT_SLASH)[indexOfSampleID] == sample.getName() } as List<File>
            List<String> listOfRunIDs = listOfFastqFilesForSample.collect { it.absolutePath.split(StringConstants.SPLIT_SLASH)[indexOfRunID] }.unique() as List<String>
            listOfRunIDs.each { String runID ->
                List<File> fastqFilesForRun = listOfFastqFilesForSample.findAll { it.absolutePath.split(StringConstants.SPLIT_SLASH)[indexOfRunID] == runID } as List<File>
                List<LaneFileGroup> bundleFiles = QCPipelineScriptFileServiceHelper.sortAndPairLaneFilesToGroupsForSampleAndRun(context, sample, libraryID, runID, fastqFilesForRun);
                laneFiles.addAll(bundleFiles);
            }

        } else {
            // Default case

            File sampleDirectory = getSampleDirectory(context, sample, libraryID);

            logger.postAlwaysInfo("Searching for lane files in directory ${sampleDirectory}")
            List<File> runsForSample = FileSystemAccessProvider.getInstance().listDirectoriesInDirectory(sampleDirectory);
            for (File run : runsForSample) {
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
                laneFiles.addAll(bundleFiles);
            }
        }
        if (laneFiles.size() == 0) {
            context.addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("There were no lane files available for sample ${sample.getName()}"));
        } else {
            determineFileDecoder(context, laneFiles);
        }
        context.setProcessingFlag(flag);
        return laneFiles;
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
                String run = split[runIndex..-2].join(StringConstants.UNDERSCORE);
                String lane = String.format("L%03d", laneID);


                BamFile bamFile = COBaseFile.constructSourceFile(BamFile, f, context, new COFileStageSettings(lane, run, sample, context.getDataSet())) as BamFile
                return bamFile;
        })
        BamFileGroup bamFileGroup = new BamFileGroup(bamFiles);
        logger.info("Found ${bamFileGroup.getFilesInGroup().size()} paired bam files for sample ${sample.getName()}");
        return bamFileGroup;
    }

    public BamFileGroup getMergedBamFilesForDataSet(ExecutionContext context) {

        ProcessingFlag flag = context.setProcessingFlag(ProcessingFlag.STORE_NOTHING);

//        File outputFolderForDataSetAndAnalysis = getOutputFolderForDataSetAndAnalysis(context.getDataSet(), context.getAnalysis());
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
