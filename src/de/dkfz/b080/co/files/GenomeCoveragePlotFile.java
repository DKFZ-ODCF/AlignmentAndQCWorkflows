package de.dkfz.b080.co.files;

/**
 *
 * @author michael
 */
public class GenomeCoveragePlotFile extends COBaseFile {

    public GenomeCoveragePlotFile(CoverageTextFile parentFile) {
        super(parentFile);
    }

    public GenomeCoveragePlotFile(BamFileGroup group) {
        super(group, group.getFilesInGroup().get(0).getFileStage());
    }

}
