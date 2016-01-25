package de.dkfz.b080.co.files;

/**
 *
 * @author michael
 */
public class FlagstatsFile extends COBaseFile {

    public FlagstatsFile(BamFile bamFile) {
        super(bamFile, bamFile.getFileStage());
    }

}
