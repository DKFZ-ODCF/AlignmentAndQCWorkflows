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
import de.dkfz.roddy.execution.io.ExecutionResult;
import de.dkfz.roddy.execution.io.ExecutionService;
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider;
import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.methods.GenericMethod;
import de.dkfz.roddy.tools.LoggerWrapper;
import groovy.transform.CompileStatic;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;

/**
 * @author michael
 */
@CompileStatic
public class LaneFileGroup extends FileGroup<LaneFile> {

    private static final LoggerWrapper logger = LoggerWrapper.getLogger(LaneFile.class.getName());

    private final String id;
    private final String run;
    private final Sample sample;
    private FastqcGroup allFastqcFiles;

    public LaneFileGroup(ExecutionContext executionContext, String id, String run, Sample sample, List<LaneFile> files) {
        super(files);
        this.id = id;
        this.run = run;
        this.sample = sample;
    }

    public Sample getSample() {
        return sample;
    }

    //This method could be a candidate for
    public FastqcGroup calcFastqcForAll() {
        return calcFastqcForAll(new LinkedHashMap<String, String>());
    }

    public FastqcGroup calcFastqcForAll(Map<String, String> parameters) {
        final boolean useSingleEndProcessing = getExecutionContext().getConfiguration().getConfigurationValues().getBoolean(AlignmentAndQCConfig.FLAG_USE_SINGLE_END_PROCESSING, false);

        LinkedList<FastqcFile> files = new LinkedList<FastqcFile>();
        for (LaneFile f : filesInGroup) {
            Tuple2<FastqcFile,TextFile> calcFastqcResult = f.calcFastqc(parameters);
            files.add(calcFastqcResult.value0);

            if (useSingleEndProcessing) {
                //Directly break after the first file if single end processing is active.
                break;
            }
        }
        allFastqcFiles = new FastqcGroup(files);
        return allFastqcFiles;
    }

    public BamFile alignAndPairSlim() {
        return alignAndPairSlim(new LinkedHashMap<>());
    }

    public BamFile alignAndPairSlim(Map<String, String> parameters) {
        ExecutionContext context = getExecutionContext();
        Configuration configuration = context.getConfiguration();
        boolean useAcceleratedHardware = configuration.getConfigurationValues().getBoolean(AlignmentAndQCConfig.FLAG_USE_ACCELERATED_HARDWARE);

        // Bad hack: Decrease the file stage level by one!
        LaneFile laneFile0 = filesInGroup.get(0).getFSDecreasedCopy();
        LaneFile laneFile1 = filesInGroup.get(1).getFSDecreasedCopy();

        String libString = configuration.getConfigurationValues().getString(COConstants.PRM_CVAL_LIBRARY);
        String sampleName = laneFile0.getSample().getName();
        String pid = context.getDataSet().getId();
        String run = laneFile0.getRunID().toString();
        String lane = laneFile0.getLaneId().toString();
        String lb = sampleName + "_" + pid + (libString.equals("addToOldLib") ? "" : "_lib2");

        String laneId0 = "RAW_SEQ_FILE_1_INDEX=" + ((COFileStageSettings) laneFile0.getFileStage()).getNumericIndex();
        String laneId1 = "RAW_SEQ_FILE_2_INDEX=" + ((COFileStageSettings) laneFile1.getFileStage()).getNumericIndex();

        final String TOOL = useAcceleratedHardware ? AlignmentAndQCConfig.TOOL_ACCELERATED_ALIGNANDPAIR_SLIM.replace(':', '_') : AlignmentAndQCConfig.TOOL_ALIGNANDPAIR_SLIM;
        BamFile bamFile = GenericMethod.callGenericTool(TOOL, laneFile0, laneFile1, parameters,
                "SAMPLE=" + sampleName, "sample=" + sampleName,
                "RUN=" + run, "run=" + run,
                "LANE=" + lane, "lane=" + lane,
                "LB=" + lb, "lb=" + lb,
                parameters,
                laneId0, laneId1);
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

}
