package de.dkfz.b080.co.files

import de.dkfz.b080.co.methods.Samtools
import de.dkfz.roddy.knowledge.files.BaseFile
import groovy.transform.CompileStatic

/**
 * Created with IntelliJ IDEA.
 * User: michael
 * Date: 28.11.12
 * Time: 10:39
 * To change this template use File | Settings | File Templates.
 */
@CompileStatic
class OnTargetCoverageTextFile extends COBaseFile {

    private OnTargetCoveragePlotFile plotFile

    OnTargetCoverageTextFile(BamFile parentFile) {
        super(new BaseFile.ConstructionHelperForManualCreation(parentFile, null, null,null,null,null,null,null))
    }

    boolean hasPlotFile() {
        return plotFile != null
    }

    OnTargetCoveragePlotFile getPlotFile() {
        return plotFile
    }

    OnTargetCoveragePlotFile createPlotFile() {
        return Samtools.createOnTargetCoveragePlot(executionContext, this)
    }
}
