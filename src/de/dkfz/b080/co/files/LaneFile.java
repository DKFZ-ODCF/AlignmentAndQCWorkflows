/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

import de.dkfz.b080.co.common.IndexID;
import de.dkfz.b080.co.common.LaneID;
import de.dkfz.b080.co.common.RunID;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.io.ExecutionResult;
import de.dkfz.roddy.execution.io.ExecutionService;
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider;
import de.dkfz.roddy.knowledge.files.ITestdataSource;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.methods.GenericMethod;
import de.dkfz.roddy.tools.LoggerWrapper;
import groovy.transform.CompileStatic;

import java.io.File;
import java.util.Map;
import java.util.logging.Level;

/**
 * @author michael
 */
@CompileStatic
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

    public Tuple2<FastqcFile,TextFile> calcFastqc(Map<String, String> parameters) {
        return GenericMethod.callGenericTool("fastqc", this, parameters);
    }

    @Override
    public boolean checkFileValidity() {
        return super.checkFileValidity();    //To change body of overridden methods use File | Settings | File Templates.
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

    public String getDecompressionString() {
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
