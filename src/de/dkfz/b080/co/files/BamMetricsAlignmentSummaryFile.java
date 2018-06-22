/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

/**
 */
public class BamMetricsAlignmentSummaryFile extends COBaseFile {
    public BamMetricsAlignmentSummaryFile(BamFile parentFile) {
        super(new ConstructionHelperForManualCreation(parentFile, null, null,null,null,null,null,null));
    }
}
