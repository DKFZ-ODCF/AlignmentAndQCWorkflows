package de.dkfz.b080.co.files;

/**
 *
 * @author michael
 */
public class AlignedSequenceFile extends COBaseFile {
    public AlignedSequenceFile(LaneFile parentFile) {
        super(parentFile, parentFile.getFileStage());
//        setAsTemporaryFile();
    }
}
