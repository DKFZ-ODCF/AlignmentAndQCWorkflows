package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.LinkedList;
import java.util.List;

/**
 * Groups various coverage text files.
 * There is no check if the coverage text files have a different type of coverage.
 */
public class CoverageTextFileGroup extends FileGroup<CoverageTextFile> {

    public CoverageTextFileGroup(List<CoverageTextFile> files) {
        super(files);
    }

    public CoverageTextFileGroup() {
        super(new LinkedList<CoverageTextFile>());
    }

    /**
     * Plots a group of coverage plot files against another group of coverage plot files.
     * The plots are create on a file vs file base. So each file in the first group against each file in the other group resulting in nxm files.
     *
     * @param group The group of files against which this group should be plotted.
     * @return A complete list of coverage plots on a file by file base.
     */
    public GenomeCoveragePlotFileGroup plotAgainst(CoverageTextFileGroup group) {
        List<GenomeCoveragePlotFile> files = new LinkedList<GenomeCoveragePlotFile>();
        for (CoverageTextFile ctf : filesInGroup) {
            GenomeCoveragePlotFileGroup coveragePlotFileGroup = ctf.plotAgainst(group);
            files.addAll(coveragePlotFileGroup.getFilesInGroup());
        }
        GenomeCoveragePlotFileGroup genomeCoveragePlotFileGroup = new GenomeCoveragePlotFileGroup(files);
        return genomeCoveragePlotFileGroup;
    }

    /**
     * Plots a group of coverage plot files. No plots against other groups is done!
     * TODO: This could lead to possible errors on rerun as maybe missing group plots are not recognized.
     */

    public GenomeCoveragePlotFileGroup plot() {
        List<GenomeCoveragePlotFile> files = new LinkedList<GenomeCoveragePlotFile>();
        for (CoverageTextFile ctf : filesInGroup) {
            files.add(ctf.plot());
        }
        return new GenomeCoveragePlotFileGroup(files);
    }
}
