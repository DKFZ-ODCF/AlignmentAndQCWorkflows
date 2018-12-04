/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

import de.dkfz.b080.co.common.COConstants;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.List;

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

}
