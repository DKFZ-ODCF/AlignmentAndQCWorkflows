/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files

import de.dkfz.b080.co.common.AlignmentAndQCConfig
import de.dkfz.b080.co.common.COProjectsRuntimeService
import de.dkfz.roddy.config.RecursiveOverridableMapContainerForConfigurationValues
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.jobs.BEJobResult
import de.dkfz.roddy.execution.jobs.Job
import de.dkfz.roddy.execution.jobs.JobManager
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.knowledge.files.FileGroup
import de.dkfz.roddy.knowledge.methods.GenericMethod
import de.dkfz.roddy.tools.LoggerWrapper

/**
 * @author michael
 */
@groovy.transform.CompileStatic
class BamFileGroup extends FileGroup<BamFile> {

    static final String MERGEANDRMDUP = "mergeAndRemoveDuplicates";
    static final String MERGEANDMORMDUP_SLIM_BIOBAMBAM = "mergeAndRemoveDuplicatesSlimBioBambam";
    static final String MERGEANDMORMDUP_SLIM_PICARD = "mergeAndRemoveDuplicatesSlimPicard";
    static final String MERGEANDMORMDUP_SLIM_SAMBAMBA = "mergeAndRemoveDuplicatesSlimSambamba";

    private static LoggerWrapper logger = LoggerWrapper.getLogger(COProjectsRuntimeService.class.getName());

    private BamFile mergedBam = null;

    BamFileGroup() {
        super(null)
    }

    BamFileGroup(List<BamFile> files) {
        super(files)
    }

    BamFile mergeSlim(Sample sample, Object... additionalMergeParameters) {
        return baseMergeAndRemoveDuplicatesSlim(sample, additionalMergeParameters + ["MERGE_BAM_ONLY": true.toString()])
    }

    BamFile mergeAndRemoveDuplicatesSlimWithLibrary(Sample sample, String library, Object... additionalMergeParameters) {
        return baseMergeAndRemoveDuplicatesSlim(sample, additionalMergeParameters + ["LIBRARY": library, "library": library])
    }

    BamFile mergeAndRemoveDuplicatesSlim(Sample sample, Object... additionalMergeParameters) {
        return baseMergeAndRemoveDuplicatesSlim(sample, additionalMergeParameters)
    }

    BamFile baseMergeAndRemoveDuplicatesSlim(Sample sample, Object... additionalMergeParameters) {
        if (mergedBam == null) {
            RecursiveOverridableMapContainerForConfigurationValues cvalues = executionContext.getConfiguration().getConfigurationValues()
            boolean useBioBamBamMarkDuplicates = cvalues.getBoolean(AlignmentAndQCConfig.FLAG_USE_BIOBAMBAM_MARK_DUPLICATES, true);
            String markDuplicatesVariant = cvalues.getString(AlignmentAndQCConfig.CVALUE_MARK_DUPLICATES_VARIANT, null);
            String toolId

            if (markDuplicatesVariant == null || markDuplicatesVariant == "") {
                logger.postSometimesInfo("${AlignmentAndQCConfig.FLAG_USE_BIOBAMBAM_MARK_DUPLICATES} is deprecated. Use ${AlignmentAndQCConfig.CVALUE_MARK_DUPLICATES_VARIANT}.")
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

            mergedBam = (BamFile) GenericMethod.callGenericTool(toolId, getFilesInGroup().get(0),
                    ([this as Object] + additionalMergeParameters.toList() << ["SAMPLE": sample.getName(), "sample": sample.getName()]) as Object[])
        }
        return mergedBam
    }

    BamFile mergeAndRemoveDuplicates() {
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
        BEJobResult jobResult = new Job(run, JobManager.getInstance().createJobName(parentFiles[0], MERGEANDRMDUP, true), MERGEANDRMDUP, parameters, parentFiles, filesToVerify).run();
        bamFile.setCreatingJobsResult(jobResult);
        bamFile.getMetricsFile().setCreatingJobsResult(jobResult);
        bamFile.getFlagstatsFile().setCreatingJobsResult(jobResult);
        if (useBioBamBamMarkDuplicates) bamFile.getIndexFile().setCreatingJobsResult(jobResult);
        return bamFile;
    }
}
