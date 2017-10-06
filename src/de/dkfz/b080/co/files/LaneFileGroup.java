package de.dkfz.b080.co.files;

import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.io.ExecutionResult;
import de.dkfz.roddy.execution.io.ExecutionService;
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider;
import de.dkfz.roddy.execution.jobs.BEJobResult;
import de.dkfz.roddy.execution.jobs.Job;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;
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
public class LaneFileGroup extends FileGroup<LaneFile> {

    private static final LoggerWrapper logger = LoggerWrapper.getLogger(LaneFile.class.getName());

    private final String id;
    private final String run;
    private final Sample sample;
    private FastqcGroup allFastqcFiles;
    private AlignedSequenceFileGroup allAlignedFiles;

    public LaneFileGroup(ExecutionContext executionContext, String id, String run, Sample sample, List<LaneFile> files) {
        super(files);
        this.id = id;
        this.run = run;
        this.sample = sample;
    }

    @Override
    public void runDefaultOperations() {
        calcFastqcForAll();
        alignAll();
    }

    public boolean hasAllAligned() {
        return allAlignedFiles != null;
    }

    public boolean hasAllFastqc() {
        return allFastqcFiles != null;
    }

    public AlignedSequenceFileGroup getAllAlignedFiles() {
        return allAlignedFiles;
    }

    public FastqcGroup getAllFastqcFiles() {
        return allFastqcFiles;
    }

    public Sample getSample() {
        return sample;
    }

    //This method could be a candidate for
    public FastqcGroup calcFastqcForAll() {
        final boolean useSingleEndProcessing = getExecutionContext().getConfiguration().getConfigurationValues().getBoolean(COConstants.FLAG_USE_SINGLE_END_PROCESSING, false);

        LinkedList<FastqcFile> files = new LinkedList<FastqcFile>();
        for (LaneFile f : filesInGroup) {
            FastqcFile calcFastqc = f.calcFastqc();
            files.add(calcFastqc);

            if (useSingleEndProcessing) {
                //Directly break after the first file if single end processing is active.
                break;
            }
        }
        allFastqcFiles = new FastqcGroup(files);
        return allFastqcFiles;
    }

    public AlignedSequenceFileGroup alignAll() {
        final boolean useSingleEndProcessing = getExecutionContext().getConfiguration().getConfigurationValues().getBoolean(COConstants.FLAG_USE_SINGLE_END_PROCESSING, false);
        LinkedList<AlignedSequenceFile> files = new LinkedList<AlignedSequenceFile>();
        int i = 0;
        for (LaneFile f : filesInGroup) {
            if (useSingleEndProcessing && i > 0) { //In single end processing, every file after the first one is handled differently.
                AlignedSequenceFile alignedSequenceFile2 = new AlignedSequenceFile(f);
                alignedSequenceFile2.setCreatingJobsResult(files.get(0).getCreatingJobsResult());
                files.add(alignedSequenceFile2);
                i++;
                continue;
            }
            AlignedSequenceFile file = f.align();
            files.add(file);
            i++;
        }
        allAlignedFiles = new AlignedSequenceFileGroup(files);
        allAlignedFiles.setParentGroup(this);
        return allAlignedFiles;
    }

    public BamFile alignAndPairSlim() {
        ExecutionContext context = getExecutionContext();
        Configuration configuration = context.getConfiguration();
        boolean useAcceleratedHardware = configuration.getConfigurationValues().getBoolean(COConstants.FLAG_USE_ACCELERATED_HARDWARE);

        // Bad hack: Decrease the file stage level by one!
        LaneFile laneFile0 = filesInGroup.get(0).getFSDecreasedCopy();
        LaneFile laneFile1 = filesInGroup.get(1).getFSDecreasedCopy();

        String libString = configuration.getConfigurationValues().getString(COConstants.PRM_CVAL_LIBRARY);
        String sampleName = laneFile0.getSample().getName();
        String pid = context.getDataSet().getId();
        String run = laneFile0.getRunID();
        String lane = laneFile0.getLaneId();
        String lb = sampleName + "_" + pid + (libString.equals("addToOldLib") ? "" : "_lib2");

        String laneId0 = "RAW_SEQ_FILE_1_INDEX=" + ((COFileStageSettings) laneFile0.getFileStage()).getNumericIndex();
        String laneId1 = "RAW_SEQ_FILE_2_INDEX=" + ((COFileStageSettings) laneFile1.getFileStage()).getNumericIndex();

        final String TOOL = useAcceleratedHardware ? COConstants.TOOL_ACCELERATED_ALIGNANDPAIR_SLIM : COConstants.TOOL_ALIGNANDPAIR_SLIM;
        BamFile bamFile = GenericMethod.callGenericTool(TOOL, laneFile0, laneFile1, "SAMPLE=" + sampleName, "RUN=" + run, "LANE=" + lane, "LB=" + lb, laneId0, laneId1);
        return bamFile;
    }

    public BamFile alignAndPair() {
        ExecutionContext run = getExecutionContext();
        LaneFile laneFile0 = filesInGroup.get(0);
        LaneFile laneFile1 = filesInGroup.get(1);

        Configuration configuration = run.getConfiguration();
        BamFile bamFile = (BamFile) BaseFile.constructManual(BamFile.class, this);
        FlagstatsFile flagstatsFile = (FlagstatsFile)BaseFile.constructManual(FlagstatsFile.class, bamFile);
        BamIndexFile bamIndexFile = (BamIndexFile) BaseFile.constructManual(BamIndexFile.class, bamFile);

        //Which info is necessary? File timestamp, maybe svn version, last changes, last file, parameters?
        String libString = configuration.getConfigurationValues().getString(COConstants.PRM_CVAL_LIBRARY);
        boolean useAdaptorTrimming = configuration.getConfigurationValues().getBoolean(COConstants.FLAG_USE_ADAPTOR_TRIMMING, false);
        boolean useAcceleratedHardware = configuration.getConfigurationValues().getBoolean(COConstants.FLAG_USE_ACCELERATED_HARDWARE);
        boolean useBioBamBamSort = configuration.getConfigurationValues().getBoolean(COConstants.FLAG_USE_BIOBAMBAM_SORT);
        boolean indexCreated = !useAcceleratedHardware;
        final String TOOL = useAcceleratedHardware ? COConstants.TOOL_ACCELERATED_ALIGNANDPAIR : COConstants.TOOL_ALIGNANDPAIR;

        String sampleName = laneFile0.getSample().getName();
        String pid = run.getDataSet().getId();

        Map<String, Object> parameters = run.getDefaultJobParameters(TOOL);
        parameters.put(COConstants.PRM_FILENAME_SORTED_BAM, bamFile.getAbsolutePath());
        parameters.put(COConstants.PRM_FILENAME_FLAGSTAT, flagstatsFile.getAbsolutePath());
        parameters.put(COConstants.PRM_FILENAME_BAM_INDEX, bamIndexFile.getAbsolutePath());
        parameters.put(COConstants.PRM_RAW_SEQ_1, laneFile0.getAbsolutePath());
        parameters.put(COConstants.PRM_RAW_SEQ_2, laneFile1.getAbsolutePath());
        parameters.put(COConstants.PRM_ID, laneFile0.getRunID() + "_" + laneFile0.getLaneId());
        parameters.put(COConstants.PRM_SM, "sample_" + sampleName + "_" + pid);
        parameters.put(COConstants.PRM_LB, sampleName + "_" + pid + (libString.equals("addToOldLib") ? "" : "_lib2"));

        if (useAdaptorTrimming) {
            parameters.put(COConstants.PRM_RAW_SEQ_FILE_1_INDEX, "" + ((COFileStageSettings) laneFile0.getFileStage()).getNumericIndex());
            parameters.put(COConstants.PRM_RAW_SEQ_FILE_2_INDEX, "" + ((COFileStageSettings) laneFile1.getFileStage()).getNumericIndex());
        }

        List<BaseFile> parentFiles = new LinkedList<>();
        parentFiles.addAll(filesInGroup);
        Job job = new Job(run, run.createJobName(parentFiles.get(0), TOOL, true), TOOL, null, parameters, parentFiles, Arrays.asList((BaseFile) bamFile, flagstatsFile));
        BEJobResult jobResult = job.run();
        flagstatsFile.setCreatingJobsResult(jobResult);
        if (indexCreated) bamIndexFile.setCreatingJobsResult(jobResult);
        if (indexCreated) bamFile.setIndexFile(bamIndexFile);
        bamFile.setCreatingJobsResult(jobResult);
        bamFile.setFlagstatsFile(flagstatsFile);
        return bamFile;
    }

    public String getId() {
        return id;
    }

    public String getRun() {
        return run;
    }

    @Override
    public String toString() {
        return "LaneFileGroup{" + "id=" + id + ", run=" + run + '}';
    }

    //TODO This should be created in a remote script.
    public void createTestDataForLaneFiles() {
        //create target directory
        ExecutionContext context = getExecutionContext();
        ExecutionService es = ExecutionService.getInstance();
        FileSystemAccessProvider fs = FileSystemAccessProvider.getInstance();
        Configuration cfg = context.getConfiguration();
        String inputBasePath = context.getInputDirectory().getAbsolutePath();
        String testdataBasePath = cfg.getConfigurationValues().get("testdataBaseDirectory").toFile(context).toString();
        String testLineCnt = cfg.getConfigurationValues().get("testdataLaneFileLineCount", "10000").toString();

        StringBuilder fullCommandList = new StringBuilder();
        for (LaneFile lf : filesInGroup) {

            String filePath = lf.getPath().getAbsolutePath();
            String fileRelativePath = filePath.substring(inputBasePath.length());

            File targetFilePath = new File(testdataBasePath + File.separator + context.getDataSet() + fileRelativePath);

            if (!fs.checkDirectory(targetFilePath.getParentFile(), context, true)) {
                throw new RuntimeException("Could not create output directory " + targetFilePath.getParentFile());
            }

            fullCommandList.append(String.format("%s %s | head -n %s | %s > %s & ", lf.getDecompressionString(), filePath, testLineCnt, lf.getRecompressionString(), targetFilePath));
        }
        fullCommandList.append(" wait");
        logger.log(Level.INFO, fullCommandList.toString());
        ExecutionResult er = es.execute(fullCommandList.toString());
        if (!er.successful) {
            logger.severe("Could not create testdata for one or more files.");
        }
    }

}
