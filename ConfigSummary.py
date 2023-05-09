# Copyright (c) 2022 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#
from collections import ChainMap
from typing import Mapping, Optional, List


def subdict_str(the_dict: Mapping[str, str], keys: List[str], default=None) \
        -> Mapping[str, Optional[str]]:
    """
    Just a helper function that allows subsetting the keys of a Dict[str,str].
    Note that key that are requested but not in the input dictionary will be returned with
    value Nane.
    """
    return {k: the_dict.get(k, default) for k in keys}


class ConfigSummary:
    """
    This encodes the interpretation of the alignment configuration and mapping to parameter sets.
    Furthermore, (and this could thus be extracted into its own class) it combines the results
    using a combination strategy.
    """

    def __init__(self):
        pass

    # Variables used by the workflow and relevant either for being reported or for deciding, which
    # values to report.

    WORKFLOW_ID = "WORKFLOW_ID"
    wes = "exomeAnalysis"
    wgbs = "bisulfiteCoreAnalysis"
    wgs = "qcAnalysis"

    PERL = "PERL_VERSION"
    PYTHON = "PYTHON_VERSION"
    R = "R_VERSION"

    BEDTOOLS = "BEDTOOLS_VERSION"
    VCFTOOLS = "VCFTOOLS_VERSION"

    SAMTOOLS = "SAMTOOLS_VERSION"
    SAMPESORT_MEMSIZE = "SAMPESORT_MEMSIZE"

    # We are not tracking the _BINARY variables, because these are derived and in some cases have
    # only temporary binding (e.g. later versions use functions to temporarily set up a binary
    # because the versions for the different steps (flagstat, dupmark, view) diverged.
    SAMBAMBA = "SAMBAMBA_VERSION"  # Only used for 'view' and 'sort'

    # Version 0.4.6 of sambamba produces the old samtools 0.1.19 compatible counts.
    SAMBAMBA_FLAGSTATS = "SAMBAMBA_FLAGSTATS_VERSION"
    SAMBAMBA_MARKDUP = "SAMBAMBA_MARKDUP_VERSION" # For duplication marking. If
    SAMBAMBA_MARKDUP_OPTIONS = "SAMBAMBA_MARKDUP_OPTS"

    JAVA = "JAVA_VERSION" # for fastqc, picard and trimmomatic
    HTSLIB = "HTSLIB_VERSION"
    PICARD = "PICARD_VERSION"
    PICARD_MARKDUP_JVM_OPTS = "PICARD_MARKDUP_JVM_OPTS"

    LIBMAUS = "LIBMAUS_VERSION" # used for biobambam
    BIOBAMBAM = "BIOBAMBAM_VERSION"

    # bamFileExists: this is an internal configuration variable!
    BWA = "BWA_VERSION"
    BWA_MEM_OPTIONS = "BWA_MEM_OPTIONS"
    BWA_MEM_THREADS = "BWA_MEM_THREADS"
    INDEX_PREFIX = "INDEX_PREFIX"

    runBwaPostAltJs = "runBwaPostAltJs"
    ALT_FILE = "ALT_FILE"
    K8_VERSION = "K8_VERSION"
    K8_BINARY = "K8_BINARY"
    bwaPostAltJsPath = "bwaPostAltJsPath"
    bwaPostAltJsMinPaRatio = "bwaPostAltJsMinPaRatio"
    bwaPostAltJsHla = "bwaPostAltJsHla"

    useExistingPairedBams = "useExistingPairedBams"    # Do not use FASTQs, but run merge
    useOnlyExistingTargetBam = "useOnlyExistingTargetBam"  # Do not even merge, just do QC
    runFastQC = "runFastQC"
    runFastQCOnly = "runFastQCOnly"
    FASTQC = "FASTQC_VERSION"  # if runFastQC = true
    useAdapterTrimming = "useAdaptorTrimming"
    # TRIMMOMATIC = "TRIMMOMATIC_VERSION"  # if useAdapterTrimming; hardcoded: 0.30
    ADAPTOR_TRIMMING_OPTIONS_0 = "ADAPTOR_TRIMMING_OPTIONS_0"
    ADAPTOR_TRIMMING_OPTIONS_1 = "ADAPTOR_TRIMMING_OPTIONS_1"
    CLIP_INDEX = "CLIP_INDEX"

    TARGET_REGIONS_FILE = "TARGET_REGIONS_FILE"
    TARGETSIZE = "TARGETSIZE"

    IS_TAGMENTATION = "IS_TAGMENTATION"
    reorderUndirectionalWGBSReadPairs = "reorderUndirectionalWGBSReadPairs"
    CYTOSINE_POSITIONS_INDEX = "CYTOSINE_POSITIONS_INDEX"
    METH_CALL_PARAMETERS = "METH_CALL_PARAMETERS"

    ON_CONVEY = "ON_CONVEY"
    # biobambam, picard, sambamba. Default: empty. If set, this option takes precedence over
    # the older useBioBamBamMarkDuplicates option.
    markDuplicatesVariant = "markDuplicatesVariant"
    # use biobambam if true, otherwise use picard
    useBioBamBamMarkDuplicates = "useBioBamBamMarkDuplicates"
    useBioBamBamSort = "useBioBamBamSort"

    # TODO # if tbiPbs then VERSION strings apply; otherwise conda env
    workflowEnvironmentScript = "workflowEnvironmentScript"
    runAlignmentOnly = "runAlignmentOnly"
    runCoveragePlots = "runCoveragePlots"
    runExomeAnalysis = "runExomeAnalysis"
    runCollectBamFileMetrics = "runCollectBamFileMetrics"

    INSERT_SIZE_LIMIT = "INSERT_SIZE_LIMIT"
    runFingerprinting = "runFingerprinting"
    fingerPrintingSitesFile = "fingerprintingSitesFile"

    runACEseqQc = "runACEseqQc"  # Turn on/off the GC correction.

    # The following are defaults for some values that may not be set and would otherwise be
    # reported as null values.
    possibly_missing_value_defaults: Mapping[str, str] = {
        runFastQCOnly: "false",
        reorderUndirectionalWGBSReadPairs: "false",
        useBioBamBamSort: "false",
        runACEseqQc: "false"
    }

    def fastqc_parameters(self,
                          plugin_versions: Mapping[str, str],
                          parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = [self.runFastQC, self.runFastQCOnly]
        if parameters.get(self.runFastQC, "false") == "true" \
                and parameters.get(self.runAlignmentOnly, "false") == "false":
            result_parameters += [self.FASTQC,
                                  self.JAVA]
        else:
            return subdict_str(parameters, result_parameters)

    def trimmomatic_parameters(self,
                               plugin_versions: Mapping[str, str],
                               parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = [self.useAdapterTrimming]
        additional_parameters = {}
        if parameters.get(self.useAdapterTrimming, "false") == "true":
            additional_parameters = {
                **additional_parameters,
                "TRIMMOMATIC_VERSION": "0.30",  # Hardcoded! Binary included!
            }
            result_parameters += [self.useAdapterTrimming,
                                  self.ADAPTOR_TRIMMING_OPTIONS_0,
                                  self.ADAPTOR_TRIMMING_OPTIONS_1,
                                  self.CLIP_INDEX,
                                  self.JAVA]
        return {**subdict_str(parameters, result_parameters),
                **additional_parameters}

    def duplication_parameters(self,
                               plugin_versions: Mapping[str, str],
                               parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = []
        picard_parameters = [self.PICARD,
                             self.PICARD_MARKDUP_JVM_OPTS,
                             self.JAVA]
        sambamba_parameters = [self.SAMBAMBA_MARKDUP,
                               self.SAMBAMBA_MARKDUP_OPTIONS]
        biobambam_parameters = [self.BIOBAMBAM,
                                self.LIBMAUS]

        # Configuration is either with markDuplicatesVariant (new, high priority) or
        # useBioBamBamMarkDuplicates (old, low priority).
        on_convey = parameters.get(self.ON_CONVEY, "false")
        mark_dup = parameters.get(self.markDuplicatesVariant, None)
        if on_convey == "true":
            result_parameters += sambamba_parameters
        else:
            if mark_dup is not None:
                result_parameters += [self.markDuplicatesVariant]
                if mark_dup == "picard":
                    result_parameters += picard_parameters
                elif mark_dup == "sambamba":
                    result_parameters += sambamba_parameters
                elif mark_dup == "biobambam":
                    result_parameters += biobambam_parameters
                else:
                    raise RuntimeError("Unknown value for " +
                                       self.markDuplicatesVariant +
                                       ": " + mark_dup)
            else:
                use_biobambam = parameters.get(self.useBioBamBamMarkDuplicates, None)
                result_parameters += [self.useBioBamBamMarkDuplicates]
                if use_biobambam is None:
                    result_parameters += picard_parameters
                else:
                    result_parameters += biobambam_parameters

        return subdict_str(parameters, result_parameters)

    def sorting_parameters(self,
                           plugin_versions: Mapping[str, str],
                           parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = []
        run_fast_qc_only = parameters.get(self.runFastQCOnly, "false")
        use_only_existing_target_bam = parameters.get(self.useOnlyExistingTargetBam, "false")
        on_convey = parameters.get(self.ON_CONVEY, "false")
        use_biobambam_sort = parameters.get(self.useBioBamBamSort, "false")
        workflow_id = parameters[self.WORKFLOW_ID]
        if run_fast_qc_only == "false" and use_only_existing_target_bam == "false":
            if workflow_id in [self.wgs, self.wes]:
                if on_convey == "true":
                    result_parameters += [self.ON_CONVEY,
                                          self.SAMBAMBA]
                elif use_biobambam_sort == "false":
                    result_parameters += [self.useBioBamBamSort,
                                          self.SAMPESORT_MEMSIZE,
                                          self.SAMTOOLS]
                elif use_biobambam_sort == "true":
                    result_parameters += [self.useBioBamBamSort,
                                          self.BIOBAMBAM]
                else:
                    raise RuntimeError(f"Unexpected value for {self.useBioBamBamSort}: " +
                                       use_biobambam_sort)
            elif workflow_id == self.wgbs:
                result_parameters += [self.useBioBamBamSort,
                                      self.SAMPESORT_MEMSIZE,
                                      self.SAMTOOLS]
        else:
            # No sorting necessary, if the BAM file exists. These will ensure the information on
            # which the output is based is tracked.
            result_parameters += [self.runFastQCOnly,
                                  self.useOnlyExistingTargetBam]

        return subdict_str(parameters, result_parameters)

    def qc_parameters(self,
                      plugin_versions: Mapping[str, str],
                      parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = [self.INSERT_SIZE_LIMIT]
        if self.SAMBAMBA_FLAGSTATS in parameters.keys():
            result_parameters += [self.SAMBAMBA_FLAGSTATS]
        else:
            result_parameters += [self.SAMTOOLS]
        return subdict_str(parameters, result_parameters)

    def alignment_parameters(self,
                             plugin_version: Mapping[str, str],
                             parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = []
        if parameters.get(self.runFastQCOnly, "false") == "false" \
                and parameters.get(self.useExistingPairedBams, "false") == "false" \
                and parameters.get(self.useOnlyExistingTargetBam, "false") == "false":
            result_parameters += [self.BWA,
                                  self.BWA_MEM_THREADS,
                                  self.BWA_MEM_OPTIONS,
                                  self.INDEX_PREFIX]
        return subdict_str(parameters, result_parameters)

    def bwa_post_alt_parameters(self,
                                plugin_version: Mapping[str, str],
                                parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = []
        if self.runBwaPostAltJs in parameters.keys() and\
                parameters.get(self.runBwaPostAltJs) == "true":
            result_parameters += [self.runBwaPostAltJs]
            if self.K8_BINARY in parameters.keys():
                result_parameters += [self.K8_BINARY]
            elif self.K8_VERSION in parameters.keys():  # explicit binary has precedence
                result_parameters += [self.K8_VERSION]
            else:
                # The default is to take the binary besides BWA executable, but this is not decided
                # by Roddy or the plugin and is not readable from the parameter file. We just take
                # a placeholder here. The unification of parameters won't happen anyway, if the
                # plugin version has changed (and with it possibly the default).
                result_parameters += ["default"]
            result_parameters += [
                self.ALT_FILE,
                self.bwaPostAltJsPath,
                self.bwaPostAltJsMinPaRatio,
                self.bwaPostAltJsHla]
        return subdict_str(parameters, result_parameters)

    def wes_parameters(self,
                       plugin_version: Mapping[str, str],
                       parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = []
        if parameters[self.WORKFLOW_ID] == self.wes:
            if parameters[self.runExomeAnalysis] != "true":
                raise RuntimeError("Inconsistent configuration: WORKFLOW_ID == 'exomeAnalysis' " +
                                   "but runExomeAnalysis != 'true'")
            result_parameters += [self.TARGET_REGIONS_FILE,
                                  self.TARGETSIZE,
                                  self.BEDTOOLS]
        return subdict_str(parameters, result_parameters)

    def wgbs_parameters(self,
                        plugin_version: Mapping[str, str],
                        parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = []
        if parameters[self.WORKFLOW_ID] == self.wgbs:
            result_parameters += [self.IS_TAGMENTATION,
                                  self.reorderUndirectionalWGBSReadPairs,
                                  self.CYTOSINE_POSITIONS_INDEX,
                                  self.METH_CALL_PARAMETERS]
        return subdict_str(parameters, result_parameters)

    def aceseqqc_parameters(self,
                            plugin_version: Mapping[str, str],
                            parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = [self.runACEseqQc]
        if parameters.get(self.runACEseqQc, "false"):
            result_parameters += [self.HTSLIB,
                                  self.VCFTOOLS]
        return subdict_str(parameters, result_parameters)

    def fingerprinting_parameters(self,
                                  plugin_version: Mapping[str, str],
                                  parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = [self.runFingerprinting]
        if parameters.get(self.runFingerprinting, "false") == "true":
            result_parameters += [self.fingerPrintingSitesFile]
        return subdict_str(parameters, result_parameters)

    def coverage_plot_parameters(self,
                                 plugin_version: Mapping[str, str],
                                 parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        result_parameters = [self.runCoveragePlots]
        if parameters.get(self.runCoveragePlots, "false"):
            result_parameters += [self.BEDTOOLS]
        return subdict_str(parameters, result_parameters)

    def summarize(self,
                  plugin_versions: Mapping[str, str],
                  parameters: Mapping[str, str]) -> \
                      Mapping[str, Mapping[str, Optional[str]]]:
        """
        This function recapitulates the logic of the workflow. Only parameters and versions are
        returned that were actually relevant for the analysis described by the parameters and
        plugin_versions. Therefore, version-specific differences may also be addressed here.

        Returns a nested dictionary structure with string values, that describe the parameters.
        """
        # Copy the parameters. ChainMap needs a MutableMapping.
        mutable_parameters = {k: v for k, v in parameters.items()}
        # Some parameters may not be set and would otherwise be reported as "null". Their defaults
        # are generally "false". We fill these values in, in case they are missing.
        parameters = ChainMap(mutable_parameters, self.possibly_missing_value_defaults)
        results: Mapping[str, Mapping[str, Optional[str]]] = {
            "workflow": {
                "id": parameters[self.WORKFLOW_ID]
            },
            "roddy": plugin_versions,
            "base": subdict_str(parameters,
                                [self.PERL,
                                 self.PYTHON,
                                 self.R,
                                 self.SAMTOOLS
                                 ]),
            "adapter_trimming": self.trimmomatic_parameters(plugin_versions, parameters),
            "fastqc": self.fastqc_parameters(plugin_versions, parameters),
            "duplication_marking": self.duplication_parameters(plugin_versions, parameters),
            "sorting": self.sorting_parameters(plugin_versions, parameters),
            "qc": self.qc_parameters(plugin_versions, parameters),
            "aceseq_qc": self.aceseqqc_parameters(plugin_versions, parameters),
            "fingerprinting": self.fingerprinting_parameters(plugin_versions, parameters),
            "alignment": self.alignment_parameters(plugin_versions, parameters),
            "bwa_post_alt": self.bwa_post_alt_parameters(plugin_versions, parameters),
            "wes": self.wes_parameters(plugin_versions, parameters),
            "wgbs": self.wgbs_parameters(plugin_versions, parameters),
        }
        return results