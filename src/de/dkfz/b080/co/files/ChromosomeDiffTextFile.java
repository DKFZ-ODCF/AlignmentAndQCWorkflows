package de.dkfz.b080.co.files;

/**
 * These files are always temporary
 * @author michael
 */
public class ChromosomeDiffTextFile extends COBaseFile {
    public ChromosomeDiffTextFile(BamFile parentFile) {
        super(parentFile);
        setAsTemporaryFile();
    }

}
