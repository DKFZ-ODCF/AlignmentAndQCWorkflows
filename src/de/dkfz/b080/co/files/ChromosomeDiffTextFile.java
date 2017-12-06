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
