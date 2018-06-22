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
public class AlignedSequenceFile extends COBaseFile {
    public AlignedSequenceFile(LaneFile parentFile) {
        super(new ConstructionHelperForGenericCreation(parentFile, null, null, null, null, null, parentFile.getFileStage(), null));
//        setAsTemporaryFile();
    }
}
