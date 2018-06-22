/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files
import de.dkfz.roddy.knowledge.files.BaseFile.ConstructionHelperForBaseFiles

@groovy.transform.CompileStatic
class ChromosomeDiffValueFile extends COBaseFile {
    ChromosomeDiffValueFile(ConstructionHelperForBaseFiles helper) {
        super(helper)
    }
}
