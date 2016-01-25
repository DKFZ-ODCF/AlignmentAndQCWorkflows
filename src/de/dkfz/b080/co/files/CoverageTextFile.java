package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.LinkedList;
import java.util.List;

/**
 * @author michael
 */
public class CoverageTextFile extends COBaseFile {

    public enum CoverageType {
        TargetEnrichment,
        RawBamTargetEnrichment,
        Default,
        ReadBins
    }

    private final CoverageType coverageType;

    public CoverageTextFile(BamFile bamFile) {
        this(bamFile, CoverageType.Default);
    }

    public CoverageTextFile(BamFile bamFile, CoverageType coverageType) {
        super(bamFile);
        this.coverageType = coverageType;
    }

    public CoverageType getCoverageType() {
        return coverageType;
    }

    /**
     * Plots a coverage file against a range of other files.
     */
    public GenomeCoveragePlotFileGroup plotAgainst(CoverageTextFileGroup coverageTextFiles) {
        List<GenomeCoveragePlotFile> files = new LinkedList<GenomeCoveragePlotFile>();
        for (CoverageTextFile ctf : coverageTextFiles.getFilesInGroup()) {
            files.add(plotAgainst(ctf));
        }
        GenomeCoveragePlotFileGroup plotFileGroup = new GenomeCoveragePlotFileGroup(files);
        return plotFileGroup;
    }

    public GenomeCoveragePlotFile plotAgainst(CoverageTextFile coverageTextFile) {
        return GenericMethod.callGenericTool("coveragePlot", (BaseFile)this, coverageTextFile);
    }

    public GenomeCoveragePlotFile plot() {
        GenomeCoveragePlotFile genomeCoveragePlotFile = GenericMethod.callGenericTool("coveragePlotSingle", (BaseFile)this);
        genomeCoveragePlotFile.overrideFilenameUsingSelectionTag("singlePlot");
        return genomeCoveragePlotFile;
    }
}
