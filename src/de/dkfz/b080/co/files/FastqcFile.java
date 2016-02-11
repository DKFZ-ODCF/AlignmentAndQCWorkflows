package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * Represents a fastqc - quality control - file
 * @author michael
 */
public class FastqcFile extends COBaseFile {

    public FastqcFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
//        super(new ConstructionHelperForManualCreation(parentFile, null, null,null,null,null,null,null));
    }
}
