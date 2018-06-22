/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

/**
 *
 * @author michael
 */
public class GenomeCoveragePlotFile extends COBaseFile {

    public GenomeCoveragePlotFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public GenomeCoveragePlotFile(CoverageTextFile parentFile) {
        super(new ConstructionHelperForManualCreation(parentFile, null, null,null,null,null,null,null));
    }

    public GenomeCoveragePlotFile(BamFileGroup group) {
        super(new ConstructionHelperForManualCreation(group, null, null,null,null,null,null,null));
    }

}
