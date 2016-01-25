package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.List;

/**
 * @author michael
 */
public class ChromosomeDiffFileGroup extends FileGroup<BaseFile> {

    private final ChromosomeDiffPlotFile plotFile;
    private final ChromosomeDiffTextFile textFile;
    private final ChromosomeDiffValueFile valueFile;

    public ChromosomeDiffFileGroup(List<BaseFile> allFiles) {
        super(allFiles);
        ChromosomeDiffTextFile tf = null;
        ChromosomeDiffPlotFile pf = null;
        ChromosomeDiffValueFile vf = null;
        for (BaseFile bf : allFiles) {
            if (bf instanceof ChromosomeDiffTextFile) {
                tf = (ChromosomeDiffTextFile) bf;
            } else if (bf instanceof ChromosomeDiffPlotFile) {
                pf = (ChromosomeDiffPlotFile) bf;
            } else if (bf instanceof ChromosomeDiffValueFile) {
                vf = (ChromosomeDiffValueFile) bf;
            }
        }
        this.plotFile = pf;
        this.textFile = tf;
        this.valueFile = vf;
    }

    public ChromosomeDiffPlotFile getPlotFile() {
        return plotFile;
    }

    public ChromosomeDiffTextFile getTextFile() {
        return textFile;
    }

    public ChromosomeDiffValueFile getValueFile() {
        return valueFile;
    }
}
