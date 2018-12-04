/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files

import de.dkfz.roddy.knowledge.files.BaseFile.ConstructionHelperForManualCreation;

/**
 * Created with IntelliJ IDEA.
 * User: michael
 * Date: 28.11.12
 * Time: 10:40
 * To change this template use File | Settings | File Templates.
 */
@groovy.transform.CompileStatic
public class OnTargetCoveragePlotFile extends COBaseFile {
    public OnTargetCoveragePlotFile(OnTargetCoverageTextFile parentFile) {
        super(new ConstructionHelperForManualCreation(parentFile, null, null,null,null,null,null,null));
    }
}
