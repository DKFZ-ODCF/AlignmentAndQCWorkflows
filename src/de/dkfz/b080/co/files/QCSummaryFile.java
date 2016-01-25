package de.dkfz.b080.co.files;

import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.b080.co.methods.Common;

import java.util.List;

/**
 * @author michael
 */
public class QCSummaryFile extends COBaseFile {

    public QCSummaryFile(BamFile bamFile) {
        super(bamFile);
    }

    public QCSummaryFile(BamFile bamFile, List<COBaseFile> parentFiles) {
        super(bamFile);
        for (COBaseFile bf : parentFiles) {
            if (bf instanceof BamFile) continue;
            this.parentFiles.add(bf);
        }
    }

    public static QCSummaryFile createFromFileList(ExecutionContext executionContext, BamFile bamFile, List<COBaseFile> files) {
        return Common.createQCSummaryFileFromList(executionContext, bamFile, files);
    }
}
