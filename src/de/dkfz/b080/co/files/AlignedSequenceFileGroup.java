/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

import de.dkfz.b080.co.common.AlignmentAndQCConfig;
import de.dkfz.b080.co.common.COConstants;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.jobs.BEJobResult;
import de.dkfz.roddy.execution.jobs.Job;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * @author michael
 */
public class AlignedSequenceFileGroup extends FileGroup<AlignedSequenceFile> {

    public AlignedSequenceFileGroup(List<AlignedSequenceFile> files) {
        super(files);
    }

    // TODO This method could also be simpler. Most of the things are configurable or could be passed via config.
    public BamFile pairAndSortSlim() {
        ExecutionContext context = getExecutionContext();
        Configuration configuration = context.getConfiguration();
      
        AlignedSequenceFile seqFile0 = filesInGroup.get(0);
        AlignedSequenceFile seqFile1 = filesInGroup.get(1);
        LaneFile laneFile0 = (LaneFile) seqFile0.getParentFiles().get(0);
        LaneFile laneFile1 = (LaneFile) seqFile1.getParentFiles().get(0);

        String libString = configuration.getConfigurationValues().getString(COConstants.PRM_CVAL_LIBRARY);
        String sampleName = laneFile0.getSample().getName();
        String pid = context.getDataSet().getId();
        String run = laneFile0.getRunID().toString();
        String lane = laneFile0.getLaneId().toString();
        String lb = sampleName + "_" + pid + (libString.equals("addToOldLib") ? "" : "_lib2");

        String laneId0 = "RAW_SEQ_FILE_1_INDEX=" + ((COFileStageSettings) laneFile0.getFileStage()).getNumericIndex();
        String laneId1 = "RAW_SEQ_FILE_2_INDEX=" + ((COFileStageSettings) laneFile1.getFileStage()).getNumericIndex();

        final String TOOL = "sampesortSlim";
        BamFile bamFile = GenericMethod.callGenericTool(TOOL, seqFile0, seqFile1, laneFile0, laneFile1, "SAMPLE=" + sampleName, "RUN=" + run, "LANE=" + lane, "LB=" + lb, laneId0, laneId1);
        return bamFile;
    }

    // TODO Rework the tool? Make the method simple. It is possible now.
    public BamFile pairAndSort() {
        ExecutionContext context = getExecutionContext();

        AlignedSequenceFile seqFile0 = filesInGroup.get(0);
        AlignedSequenceFile seqFile1 = filesInGroup.get(1);
        LaneFile laneFile0 = (LaneFile) seqFile0.getParentFiles().get(0);
        LaneFile laneFile1 = (LaneFile) seqFile1.getParentFiles().get(0);

        LaneFile parentFile = (LaneFile) filesInGroup.get(0).getParentFiles().get(0);
        Configuration configuration = context.getConfiguration();

        // Can be configured
        BamFile bamFile = (BamFile)BaseFile.constructManual(BamFile.class, this);
        FlagstatsFile flagstatsFile = (FlagstatsFile)BaseFile.constructManual(FlagstatsFile.class, bamFile);
        BamIndexFile indexFile = (BamIndexFile) BaseFile.constructManual(BamIndexFile.class, bamFile);

        // Which info is necessary? File timestamp, maybe svn version, last changes, last file, parameters?
        String libString = configuration.getConfigurationValues().get(COConstants.PRM_CVAL_LIBRARY).toString();
        boolean useAdaptorTrimming = configuration.getConfigurationValues().getBoolean(AlignmentAndQCConfig.FLAG_USE_ADAPTOR_TRIMMING, false);

        // All parameters via config
        Map<String, Object> parameters = context.getDefaultJobParameters( AlignmentAndQCConfig.TOOL_SAMPESORT);
        parameters.put(COConstants.PRM_FILENAME_SORTED_BAM, bamFile.getAbsolutePath());
        parameters.put(COConstants.PRM_FILENAME_SEQ_1, seqFile0.getAbsolutePath());
        parameters.put(COConstants.PRM_FILENAME_SEQ_2, seqFile1.getAbsolutePath());
        parameters.put(COConstants.PRM_FILENAME_FLAGSTAT, flagstatsFile.getAbsolutePath());
        parameters.put(COConstants.PRM_RAW_SEQ_1, laneFile0.getAbsolutePath());
        parameters.put(COConstants.PRM_RAW_SEQ_2, laneFile1.getAbsolutePath());

        // Could be moved to the config / scripts
        parameters.put("ID", parentFile.getRunID() + "_" + parentFile.getLaneId());
        parameters.put("SM", "sample_" + parentFile.getSample().getName() + "_" + context.getDataSet());
        parameters.put("LB", parentFile.getSample().getName() + "_" + context.getDataSet() + (libString.equals("addToOldLib") ? "" : "_lib2"));

        // Can be put in always.
        if(useAdaptorTrimming) {
            parameters.put(COConstants.PRM_RAW_SEQ_FILE_1_INDEX, "" + ((COFileStageSettings)seqFile0.getFileStage()).getNumericIndex());
            parameters.put(COConstants.PRM_RAW_SEQ_FILE_2_INDEX, "" + ((COFileStageSettings)seqFile1.getFileStage()).getNumericIndex());
        }

        List<BaseFile> parentFiles = new LinkedList<BaseFile>();
        parentFiles.addAll(filesInGroup);
        Job job = new Job(context, context.createJobName(parentFiles.get(0), AlignmentAndQCConfig.TOOL_SAMPESORT, true), AlignmentAndQCConfig.TOOL_SAMPESORT, null, parameters, parentFiles, Arrays.asList((BaseFile)bamFile, indexFile, flagstatsFile));
        BEJobResult jobResult = job.run();

        flagstatsFile.setCreatingJobsResult(jobResult);
        indexFile.setCreatingJobsResult(jobResult);
        bamFile.setCreatingJobsResult(jobResult);
        bamFile.setFlagstatsFile(flagstatsFile);
        bamFile.setIndexFile(indexFile);
        return bamFile;

    }

}
