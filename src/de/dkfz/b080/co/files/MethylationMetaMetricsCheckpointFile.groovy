/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files
import de.dkfz.roddy.knowledge.files.*

@groovy.transform.CompileStatic
class MethylationMetaMetricsCheckpointFile extends COBaseFile {
    public MethylationMetaMetricsCheckpointFile(BaseFile.ConstructionHelperForBaseFiles helper) {
        super(helper)
    }
}
