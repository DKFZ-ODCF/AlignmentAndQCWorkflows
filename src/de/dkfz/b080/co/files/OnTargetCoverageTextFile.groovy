/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files

import de.dkfz.b080.co.methods.Samtools
import de.dkfz.roddy.knowledge.files.BaseFile.ConstructionHelperForManualCreation;

/**
 * Created with IntelliJ IDEA.
 * User: michael
 * Date: 28.11.12
 * Time: 10:39
 * To change this template use File | Settings | File Templates.
 */
@groovy.transform.CompileStatic
public class OnTargetCoverageTextFile extends COBaseFile {

    private de.dkfz.b080.co.files.OnTargetCoveragePlotFile plotFile;

    public OnTargetCoverageTextFile(BamFile parentFile) {
        super(new ConstructionHelperForManualCreation(parentFile, null, null,null,null,null,null,null));
    }
//    OnTargetCoverageTextFile(File path, ExecutionContext executionContext, JobResult creatingJobsResult, List<BaseFile> parentFiles, FileStageSettings settings) {
//        super(path, executionContext, creatingJobsResult, parentFiles, settings)
//    }

    public boolean hasPlotFile() {
        return plotFile != null;
    }

    public de.dkfz.b080.co.files.OnTargetCoveragePlotFile getPlotFile() {
        return plotFile
    }

    public de.dkfz.b080.co.files.OnTargetCoveragePlotFile createPlotFile() {
        return Samtools.createOnTargetCoveragePlot(executionContext, this);
    }
}
