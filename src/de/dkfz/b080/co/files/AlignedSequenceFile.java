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
