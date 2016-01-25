package de.dkfz.b080.co.methods

import de.dkfz.b080.co.files.BamFile
import de.dkfz.b080.co.files.BamIndexFile
import de.dkfz.b080.co.files.OnTargetCoveragePlotFile
import de.dkfz.b080.co.files.OnTargetCoverageTextFile
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.jobs.Job
import de.dkfz.roddy.execution.jobs.JobResult
import de.dkfz.roddy.execution.jobs.ScriptCallingMethod
import de.dkfz.roddy.execution.jobs.StaticScriptProviderClass
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.knowledge.methods.GenericMethod
import de.dkfz.roddy.tools.LoggerWrapper

import java.util.logging.Logger

/**
 *
 * @author michael
 */
@groovy.transform.CompileStatic
@StaticScriptProviderClass
class Samtools {

    private static final LoggerWrapper logger = LoggerWrapper.getLogger(Samtools.class.name);

    private static final String ONTARGETCOVERAGEPLOT = "onTargetCoveragePlotter";
    private static final String ONTARGETCOVERAGEPLOTTER_BINARY = "onTargetCoveragePlotterBinary";

    @ScriptCallingMethod
    public static OnTargetCoveragePlotFile createOnTargetCoveragePlot(ExecutionContext run, OnTargetCoverageTextFile textFile) {
        OnTargetCoveragePlotFile onTargetCoveragePlotFile = new OnTargetCoveragePlotFile(textFile);
        File plotfile = onTargetCoveragePlotFile.getPath();

        Map<String, Object> parameters = new HashMap<String, Object>();
        parameters.put("TARGETS_WITH_COVERAGE_TEXT", textFile.path.absolutePath);
        parameters.put("TARGETS_PLOT", plotfile.absolutePath);
        parameters.put("FILENAME_PREFIX", textFile.fileStage.getIDString());
        parameters.put(ONTARGETCOVERAGEPLOTTER_BINARY, run.getConfiguration().getProcessingToolPath(run, ONTARGETCOVERAGEPLOTTER_BINARY).absolutePath);

        JobResult jobResult = new Job(run, run.createJobName(textFile, ONTARGETCOVERAGEPLOT), ONTARGETCOVERAGEPLOT, parameters, [(BaseFile) textFile]).run();
        onTargetCoveragePlotFile.setCreatingJobsResult(jobResult);
        return onTargetCoveragePlotFile;
    }

}
