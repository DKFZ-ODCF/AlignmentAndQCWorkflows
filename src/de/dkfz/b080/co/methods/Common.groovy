package de.dkfz.b080.co.methods
import de.dkfz.b080.co.files.*
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.jobs.Job
import de.dkfz.roddy.execution.jobs.JobResult
import de.dkfz.roddy.execution.jobs.ScriptCallingMethod
import de.dkfz.roddy.execution.jobs.StaticScriptProviderClass
import de.dkfz.roddy.knowledge.files.BaseFile
/**
 *
 * @author michael
 */
@groovy.transform.CompileStatic
@StaticScriptProviderClass
class Common {

//    public static final String FASTQC = "fastqc";
//    public static final String INSERTSIZES = "insertSizes";
    public static final String CHROMOSOMEDIFF = "chromosomeDiff";
    public static final String GENOMECOVERAGE = "genomeCoverage";
    public static final String QCSUMMARY = "qcSummary";

    public static final String CONFIG_FILE = "CONFIG_FILE";
    public static final String PID = "DataSet";
    public static final String TOOLS_DIR = "TOOLS_DIR";
    public static final String ANALYSIS_DIR = "ANALYSIS_DIR";
//
//    public static String createJobName(BaseFile bf, String toolName, boolean reduceLevel = false, List<BaseFile> inputFilesForSizeCalculation = []) {
//        return Roddy.commandFactory.createJobName(bf, toolName, reduceLevel);
//    }


//    public static void determineSequencerID(RunningProcess context, List<LaneFileGroup> laneFileGroupList) {
//        ExecutionService es = context.getExecutionService();
//        File sequencerIDTool = context.getConfiguration().getProcessingToolPath(context, "sequencerDetection");
//        List<LaneFile> allLaneFiles = [];
//        StringBuilder command = new StringBuilder();
//        for (LaneFileGroup lfg : laneFileGroupList) {
//            List<LaneFile> lflist = lfg.getFilesInGroup();
//            for (LaneFile lf : lflist) {
//                allLaneFiles << lf;
//                //TODO Compressor options!
//                command << "&& " << "${lf.decompressionString} ${lf.path.absolutePath} | ${sequencerIDTool.absolutePath} ";
//            }
//        }
//        ExecutionResult er = es.execute(command[2..-1].toString());
//        for (int i = 0; i < er.resultLines.size(); i++) {
//            allLaneFiles[i].sequencerID = er.resultLines[i];
//        }
//    }

//    public static FastqcFile fastqc(ExecutionContext context, LaneFile laneFile) {
//        FastqcFile fastqcFile = new FastqcFile(laneFile);
//        File filePath = fastqcFile.getPath();
//
//        Map<String, Object> parameters = new HashMap<String, Object>();
//        parameters.putAll([
//                PRM_RAW_SEQ: laneFile.path.absolutePath,
//                "FILENAME_FASTQC": filePath.getAbsolutePath(),
//        ]);
//
//        JobResult jobResult = new Job(context, createJobName(laneFile, FASTQC), [filePath], FASTQC, parameters).context();
//        fastqcFile.setCreatingJobsResult(jobResult);
//        return fastqcFile;
//    }

//    public static InsertSizesFileGroup determineInsertSizesForBamFile(ExecutionContext context, BamFile bamFile) {
//        if (!bamFile.hasIndex()) bamFile.indexFile();
//        InsertSizesTextFile tFile = new InsertSizesTextFile(bamFile);
//        InsertSizesPlotFile pFile = new InsertSizesPlotFile(bamFile);
//        File filePathD = tFile.path;
//        File filePathP = pFile.path;
//
//        //Find first position of .bam. Take everything from the left of this position.
//        String bamFilename = bamFile.path.absolutePath;
////        String bamFileShortname = bamFile.path.name;
//
////        String fileInfo = bamFileShortname[0..bamFileShortname.lastIndexOf("_") - 1];
//
//        Map<String, String> parameters = [
//                "FILENAME": bamFilename,
//                "FILENAMED": filePathD.absolutePath,
//                "FILENAMEP": filePathP.absolutePath
//        ];
//        List<BaseFile> pFiles = [(BaseFile) bamFile];
//        JobResult jobResult = new Job(context, createJobName(pFiles[0], INSERTSIZES), INSERTSIZES, parameters, pFiles).context();
//        tFile.setCreatingJobsResult(jobResult);
//        pFile.setCreatingJobsResult(jobResult);
//
//        InsertSizesFileGroup fGroup = new InsertSizesFileGroup(tFile, pFile, context);
//
//        return fGroup;
//    }

    @ScriptCallingMethod
    public static ChromosomeDiffFileGroup differentiateChromosomesForBamFile(ExecutionContext run, BamFile bamFile) {
        if (!bamFile.hasIndex()) bamFile.index();

        ChromosomeDiffTextFile tFile = new ChromosomeDiffTextFile(bamFile);
        ChromosomeDiffPlotFile pFile = new ChromosomeDiffPlotFile(bamFile);

        File filePathD = tFile.path;
        File filePathP = pFile.path;

        String bamFilename = bamFile.path.absolutePath;

        Map<String, Object> parameters = new HashMap<String, Object>();
        parameters.putAll([
                "FILENAME": bamFilename,
                "FILENAMED": filePathD.absolutePath,
                "FILENAMEP": filePathP.absolutePath
        ]);
        List<BaseFile> pFiles = [(BaseFile) bamFile.getIndexFile()];

        JobResult jobResult = new Job(run, run.createJobName(pFiles[0], CHROMOSOMEDIFF), CHROMOSOMEDIFF, parameters, pFiles).run();
        tFile.setCreatingJobsResult(jobResult);
        pFile.setCreatingJobsResult(jobResult);
        ChromosomeDiffFileGroup fGroup = new ChromosomeDiffFileGroup([tFile, pFile] as List<BaseFile>);
        return fGroup;
    }

//
//    @ScriptCallingMethod
//    public static CoverageTextFile calculateBamCoverage(ExecutionContext run, BamFile bamFile, CoverageTextFile.CoverageType coverageType = CoverageTextFile.CoverageType.Default) {
//        if (!bamFile.hasIndex()) bamFile.indexFile();
//        CoverageTextFile cTextFile = new CoverageTextFile(bamFile, coverageType);
//
//        Map<String, Object> parameters = new HashMap<String, Object>();
//        parameters.putAll([
//                "FILENAME_COVERAGE": cTextFile.absolutePath,
//                "FILENAME": bamFile.path.absolutePath,
//                "DIR_COVERAGE": cTextFile.containingFolder,
//                "COVERAGE_TYPE": coverageType.toString(),
//        ]);
//        List<BaseFile> pFiles = [(BaseFile) bamFile.getIndexFile()];
//        JobResult jobResult = new Job(run, run.createJobName(pFiles[0], GENOMECOVERAGE), GENOMECOVERAGE, parameters, pFiles, [(BaseFile)cTextFile]).run();
//        cTextFile.setCreatingJobsResult(jobResult);
//        return cTextFile;
//    }

    private static LaneFile recursivelySearchLaneFile(List<BaseFile> files) {
        LaneFile lf = null;
        for (BaseFile bf : files) {
            if (bf instanceof LaneFile) {
                lf = (LaneFile) bf;
            } else {
                lf = recursivelySearchLaneFile(bf.getParentFiles());
                if (lf != null) {
                    break;
                }
            }
        }
        return lf;
    }

    @ScriptCallingMethod
    public static QCSummaryFile createQCSummaryFileFromList(ExecutionContext run, BamFile bamFile, List<COBaseFile> files) {
        QCSummaryFile qcSummaryFile = new QCSummaryFile(bamFile, files);//filePath, context, jobResult, files, files[0].getFileStage());

        LaneFile laneFile = recursivelySearchLaneFile([(BaseFile) bamFile]);
        COFileStageSettings bamFileFileStage = (COFileStageSettings) bamFile.getFileStage()

        String sample;
        if(laneFile != null) {
            sample = laneFile.getSample().getName();
        } else {
            sample = bamFileFileStage.getSample().getName();
        }
        String runId = "all_merged";
        String lane = "all_merged";

        if (bamFileFileStage.stage.isMoreDetailedOrEqualTo(COFileStage.RUN)) {
            runId = bamFileFileStage.runID;
            if (bamFileFileStage.stage.isMoreDetailedOrEqualTo(COFileStage.LANE))
                lane = bamFileFileStage.laneId;

        }
        def temp = run.getDefaultJobParameters(QCSUMMARY);
        Map<String, Object> parameters = (Map<String, Object>)temp;
        parameters.putAll([
                "SAMPLE": sample,
                "RUN": runId,
                "LANE": lane,
                "FILENAME_QCSUM": qcSummaryFile.absolutePath
        ]);
        files.each {
            it ->
                BaseFile bf = (BaseFile) it;
                if (bf instanceof FlagstatsFile) {
                    parameters["FILENAME_FLAGSTAT"] = bf.path.absolutePath;
                } else if (bf instanceof ChromosomeDiffValueFile) {
                    parameters["FILENAME_DIFFCHROM"] = bf.path.absolutePath;
                } else if (bf instanceof InsertSizesValueFile) {
                    parameters["FILENAME_ISIZE"] = bf.path.absolutePath;
                } else if (bf instanceof BamMetricsFile) {
                    parameters["FILENAME_METRICS"] = bf.path.absolutePath;
//                } else if (bf instanceof BamFile && ((BamFile) bf).isTargetExtractedBamFile()) {
                } else if (bf instanceof CoverageTextFile && bamFile.getTargetCoverageTextFile() == bf) {
                    parameters["FILENAME_TARGETCAPTURE"] = bf.path.absolutePath;
                } else if (bf instanceof CoverageTextFile) {
                    parameters["FILENAME_COVERAGE"] = bf.path.absolutePath;
                }
        }

        JobResult jobResult = new Job(run, run.createJobName(files[0], QCSUMMARY), QCSUMMARY, parameters, new LinkedList<BaseFile>(files), [(COBaseFile) qcSummaryFile] as List<BaseFile>).run();
        qcSummaryFile.setCreatingJobsResult(jobResult);
        return qcSummaryFile;
    }


    public static boolean isFileValid(BaseFile basefile) {

    }
}
