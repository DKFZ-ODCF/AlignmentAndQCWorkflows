package de.dkfz.b080.co.files

import de.dkfz.b080.co.common.COConfig
import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.roddy.config.OnScriptParameterFilenamePattern
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.jobs.JobManager
import de.dkfz.roddy.execution.jobs.Job
import de.dkfz.roddy.execution.jobs.JobResult
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.knowledge.files.FileGroup
import de.dkfz.roddy.knowledge.methods.GenericMethod
import de.dkfz.roddy.tools.LoggerWrapper

/**
 * @author michael
 */
@groovy.transform.CompileStatic
public class BamFileGroup extends FileGroup<BamFile> {

    public static final String MERGEANDRMDUP = "mergeAndRemoveDuplicates";
    public static final String MERGEANDMORMDUP_SLIM_BIOBAMBAM = "mergeAndRemoveDuplicatesSlimBioBambam";
    public static final String MERGEANDMORMDUP_SLIM_PICARD = "mergeAndRemoveDuplicatesSlimPicard";
    public static final String MERGEANDMORMDUP_SLIM_SAMBAMBA = "mergeAndRemoveDuplicatesSlimSambamba";
    public static final String COVERAGEPLOT = "coveragePlots";
    public static final String SNPCOMP = "snpComparison";

    private static LoggerWrapper logger = LoggerWrapper.getLogger(COProjectsRuntimeService.class.getName());

    private BamFile mergedBam = null;

    public BamFileGroup() {
        super(null);
    }

    public BamFileGroup(List<BamFile> files) {
        super(files);
    }

    public BamFile mergeSlim(Sample sample) {
        return baseMergeAndRemoveDuplicatesSlim(sample, "MERGE_BAM_ONLY=true")
    }

    public BamFile mergeAndRemoveDuplicatesSlimWithLibrary(Sample sample, String library) {
        return baseMergeAndRemoveDuplicatesSlim(sample, "LIBRARY=${library}")
    }

    public BamFile mergeAndRemoveDuplicatesSlim(Sample sample, String... additionalMergeParameters) {
        return baseMergeAndRemoveDuplicatesSlim(sample, additionalMergeParameters);
    }

    public BamFile baseMergeAndRemoveDuplicatesSlim(Sample sample, String... additionalMergeParameters) {
        if (mergedBam == null) {
            RecursiveOverridableMapContainerForConfigurationValues cvalues = executionContext.getConfiguration().getConfigurationValues()
            boolean useBioBamBamMarkDuplicates = cvalues.getBoolean(COConstants.FLAG_USE_BIOBAMBAM_MARK_DUPLICATES, true);
            String markDuplicatesVariant = cvalues.getString(COConstants.CVALUE_MARK_DUPLICATES_VARIANT, null);
            String toolId

            if (markDuplicatesVariant == null || markDuplicatesVariant == "") {
                logger.postSometimesInfo("${COConstants.FLAG_USE_BIOBAMBAM_MARK_DUPLICATES} is deprecated. Use ${COConstants.CVALUE_MARK_DUPLICATES_VARIANT}.")
                toolId = useBioBamBamMarkDuplicates ? MERGEANDMORMDUP_SLIM_BIOBAMBAM : MERGEANDMORMDUP_SLIM_PICARD
            } else {
                switch (markDuplicatesVariant.toLowerCase()) {
                    case "biobambam":
                        toolId = MERGEANDMORMDUP_SLIM_BIOBAMBAM
                        break
                    case "picard":
                        toolId = MERGEANDMORMDUP_SLIM_PICARD
                        break
                    case "sambamba":
                        toolId = MERGEANDMORMDUP_SLIM_SAMBAMBA
                        break
                    default:
                        throw new RuntimeException("markDuplicatesVariant=${markDuplicatesVariant} is not supported")
                }
            }

            Object[] parameterList = ([this, "SAMPLE=${sample.getName()}"] as ArrayList<Object>) + (additionalMergeParameters as ArrayList<Object>) as Object[]
            mergedBam = (BamFile) GenericMethod.callGenericTool(toolId, getFilesInGroup().get(0), parameterList)
        }
        return mergedBam;
    }

    public BamFile mergeAndRemoveDuplicates() {
        ExecutionContext run = executionContext;
        boolean useBioBamBamMarkDuplicates = executionContext.getConfiguration().getConfigurationValues().getBoolean("useBioBamBamMarkDuplicates", true);

        List<BaseFile> parentFiles = new LinkedList<BaseFile>();
        getFilesInGroup().each { BamFile it -> parentFiles.add((BaseFile) it); }

        BamFile bamFile = (BamFile)BaseFile.constructManual(BamFile.class, this);
        bamFile.setFlagstatsFile((FlagstatsFile)BaseFile.constructManual(FlagstatsFile.class, bamFile));
        bamFile.setMetricsFile((BamMetricsFile)BaseFile.constructManual(BamMetricsFile.class, bamFile));
        if (useBioBamBamMarkDuplicates) bamFile.setIndexFile((BamIndexFile)BaseFile.constructManual(BamIndexFile, bamFile));
        String parentFilePaths = parentFiles.collect { BaseFile pf -> pf.path.absolutePath }.join(":");

        Map<String, Object> parameters = run.getDefaultJobParameters(MERGEANDRMDUP);
        parameters.putAll([
                "FILENAME"          : bamFile.absolutePath,
                "FILENAME_FLAGSTATS": bamFile.getFlagstatsFile().absolutePath,
                "FILENAME_METRICS"  : bamFile.getMetricsFile().absolutePath,
                "INPUT_FILES"       : parentFilePaths
        ]);

        def filesToVerify = [bamFile, bamFile.flagstatsFile, bamFile.metricsFile] as List<BaseFile>
        if (bamFile.hasIndex())
            filesToVerify << bamFile.indexFile;
        JobResult jobResult = new Job(run, JobManager.getInstance().createJobName(parentFiles[0], MERGEANDRMDUP, true), MERGEANDRMDUP, parameters, parentFiles, filesToVerify).run();
        bamFile.setCreatingJobsResult(jobResult);
        bamFile.getMetricsFile().setCreatingJobsResult(jobResult);
        bamFile.getFlagstatsFile().setCreatingJobsResult(jobResult);
        if (useBioBamBamMarkDuplicates) bamFile.getIndexFile().setCreatingJobsResult(jobResult);
        return bamFile;
    }
}
