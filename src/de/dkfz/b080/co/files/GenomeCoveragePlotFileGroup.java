/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.List;

/**
 */
public class GenomeCoveragePlotFileGroup extends FileGroup<GenomeCoveragePlotFile> {

    public GenomeCoveragePlotFileGroup(List<GenomeCoveragePlotFile> files) {
        super(files);
    }
}
