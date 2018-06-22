/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * These files are always temporary
 * @author michael
 */
public class ChromosomeDiffTextFile extends COBaseFile {
    public ChromosomeDiffTextFile(BaseFile.ConstructionHelperForBaseFiles helper) {
        super(helper);
        setAsTemporaryFile();
    }
}
