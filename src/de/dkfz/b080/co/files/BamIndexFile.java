package de.dkfz.b080.co.files;

import de.dkfz.roddy.execution.jobs.JobResult;

/**
 *
 * @author michael
 */
public class BamIndexFile extends COBaseFile {

    public BamIndexFile(BamFile bamFile) {
        super(bamFile, bamFile.getFileStage());
    }

}
