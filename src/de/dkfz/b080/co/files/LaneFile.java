package de.dkfz.b080.co.files;

import de.dkfz.b080.co.common.IndexID;
import de.dkfz.b080.co.common.LaneID;
import de.dkfz.b080.co.common.RunID;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.io.ExecutionResult;
import de.dkfz.roddy.execution.io.ExecutionService;
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider;
import de.dkfz.roddy.execution.jobs.Job;
import de.dkfz.roddy.execution.jobs.JobResult;
import de.dkfz.roddy.knowledge.files.*;
import de.dkfz.roddy.knowledge.methods.GenericMethod;
import de.dkfz.roddy.tools.LoggerWrapper;

import java.io.File;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;

/**
 * @author michael
 */
@groovy.transform.CompileStatic
public class LaneFile extends COBaseFile implements ITestdataSource {
    private static final LoggerWrapper logger = LoggerWrapper.getLogger(LaneFile.class.getName());

    private String sequencerID = null;
    private String decompressionString = null;
    private String recompressionString = null;
    private FastqcFile fastqcFile;
    private AlignedSequenceFile alignedSequenceFile;

    public LaneFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    /** Copy constructor **/
    public LaneFile(LaneFile parent, ExecutionContext newContext) {
        super(parent);
        this.sequencerID = parent.sequencerID;
        this.decompressionString = parent.decompressionString;
        this.recompressionString = parent.recompressionString;
        this.fastqcFile = parent.fastqcFile;
        this.alignedSequenceFile = parent.alignedSequenceFile;
        this.setExecutionContext(newContext);
    }

    /**
     * Hack method! Allow to decrease a file stage
     */
    public LaneFile getFSDecreasedCopy() {
        LaneFile lf = new LaneFile(this, this.getExecutionContext());
        lf.fileStageSettings = fileStageSettings.decreaseLevel();
        return lf;
    }

    @Override
    public void runDefaultOperations() {
        calcFastqc();
        align();
    }

    public FastqcFile calcFastqc() {
        fastqcFile = GenericMethod.callGenericTool("fastqc", this);//Common.fastqc(executionContext, this);
        return fastqcFile;
    }

    @Override
    public boolean checkFileValidity() {
        return super.checkFileValidity();    //To change body of overridden methods use File | Settings | File Templates.
    }

    public AlignedSequenceFile align() {
        ExecutionContext context = getExecutionContext();
        AlignedSequenceFile alignedSequenceFile = new AlignedSequenceFile(this);
        Configuration configuration = context.getConfiguration();
        boolean useAcceleratedHardware = configuration.getConfigurationValues().getBoolean(COConstants.FLAG_USE_ACCELERATED_HARDWARE);
        boolean useAdaptorTrimming = configuration.getConfigurationValues().getBoolean(COConstants.FLAG_USE_ADAPTOR_TRIMMING, false);

        List<BaseFile> pFiles = new LinkedList<>();

        final String TOOL = useAcceleratedHardware ? COConstants.TOOL_ACCELERATED_ALIGNMENT : COConstants.TOOL_ALIGNMENT;

        Map<String, Object> parameters = context.getDefaultJobParameters(TOOL);
        parameters.put(COConstants.PRM_RAW_SEQ, getPath().getAbsolutePath());
        parameters.put(COConstants.PRM_RAW_SEQ_JOBJ, this);
        parameters.put(COConstants.PRM_FILENAME_ALIGNMENT, alignedSequenceFile.getAbsolutePath());

        if (useAdaptorTrimming) {

            LaneFile sisterFile = null;
            for (FileGroup fg : (List<FileGroup>) getFileGroups()) {
                if (fg instanceof LaneFileGroup) {
                    for (LaneFile laneFile : ((LaneFileGroup) fg).getFilesInGroup()) {
                        if (laneFile == this)
                            continue;
                        sisterFile = laneFile;
                        break;
                    }
                    break;
                }
            }

            COFileStageSettings fs = (COFileStageSettings) this.getFileStage();
            COFileStageSettings fss = (COFileStageSettings) sisterFile.getFileStage();

            parameters.put(COConstants.PRM_RAW_SEQ_FILE_1_INDEX, "" + fs.getNumericIndex());
            parameters.put(COConstants.PRM_RAW_SEQ_2, sisterFile.getPath().getAbsolutePath());
            parameters.put(COConstants.PRM_RAW_SEQ_FILE_2_INDEX, "" + fss.getNumericIndex());
        }

        Job job = new Job(context, context.createJobName(this, TOOL), TOOL, parameters, pFiles, Arrays.asList((BaseFile) alignedSequenceFile));

        JobResult jobResult = job.run();
        alignedSequenceFile.setCreatingJobsResult(jobResult);

        this.alignedSequenceFile = alignedSequenceFile;
        return alignedSequenceFile;
    }

    public boolean isAligned() {
        return alignedSequenceFile != null;
    }

    public boolean hasQualityControlFile() {
        return fastqcFile != null;
    }

    public FastqcFile getFastqcFile() {
        return fastqcFile;
    }

    public AlignedSequenceFile getAlignedSequenceFile() {
        return alignedSequenceFile;
    }

    public LaneID getLaneId() {
        return ((COFileStageSettings) fileStageSettings).getLaneId();
    }

    public IndexID getIndex() {
        return ((COFileStageSettings) fileStageSettings).getIndex();
    }

    public RunID getRunID() {
        return ((COFileStageSettings) fileStageSettings).getRunID();
    }

    public Sample getSample() {
        return ((COFileStageSettings) fileStageSettings).getSample();
    }
//
//    public String getSequencerID() {
//        //TODO Read out sequencer ID on demand? or instantly?
//        if (sequencerID == null) {
//            //Read it out and store it.
//            File tool = getRunningProcess().getConfiguration().getProcessingToolPath(getRunningProcess(), "sequencerDetection");
//            ExecutionResult er = getRunningProcess().getExecutionService().execute(this.decompressionString + "  " + this.path.getAbsolutePath() + " | " + tool.getAbsolutePath() + "");
//            if (er.successful)
//                sequencerID = er.firstLine;
//        }
//        return sequencerID;
//    }

    public void setSequencerID(String id) {
        sequencerID = id;
    }

    public String getDecompressionString() {
//        if(decompressionString == null) {
//            decompressionString = executionContext.getRuntimeService().getCompressorForLaneFile(executionContext, this);
//        }
        return decompressionString;
    }

    public String getRecompressionString() {
        return recompressionString;
    }

    public void setDecompressionString(String decompressionString) {
        this.decompressionString = decompressionString;
    }

    public void setRecompressionString(String recompressionString) {
        this.recompressionString = recompressionString;
    }

    @Override
    public String toString() {
        return "LaneFile{" + "sequencerID=" + sequencerID + ", fastqcFile=" + fastqcFile + ", alignedSequenceFile=" + alignedSequenceFile + '}';
    }

    @Override
    public boolean createTestData() {
        //create target directory
        ExecutionContext context = getExecutionContext();
        Configuration cfg = context.getConfiguration();
        String inputBasePath = context.getInputDirectory().getAbsolutePath();
        String testdataBasePath = cfg.getConfigurationValues().get("testdataBaseDirectory").toFile(context).toString();
        String filePath = path.getAbsolutePath();
        String testLineCnt = cfg.getConfigurationValues().get("testdataLaneFileLineCount", "10000").toString();

        String fileRelativePath = filePath.substring(inputBasePath.length());
        File targetFilePath = new File(testdataBasePath + File.separator + context.getDataSet() + fileRelativePath);
        ExecutionService es = ExecutionService.getInstance();

        if (!FileSystemAccessProvider.getInstance().checkDirectory(targetFilePath.getParentFile(), context, true)) {
            throw new RuntimeException("Could not create output directory " + targetFilePath.getParentFile());
        }
//        String zipTool = cfg.getConfigurationValue("ZIPTOOL").toString();
//        String zipToolOptions = cfg.getConfigurationValue("ZIPTOOL_OPTIONS").toString();
        String cmd = String.format(getDecompressionString() + " %s | head -n %s | " + recompressionString + " > %s", filePath, testLineCnt, targetFilePath);
        logger.log(Level.INFO, cmd);
        ExecutionResult er = es.execute(cmd);
        if (!er.successful) {
            logger.severe("Could not create testdata for file " + targetFilePath);
        }

        return targetFilePath != null;
    }
}
