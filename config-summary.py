#!/bin/env python3
#
# Copyright (c) 2021 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#
# $ config-summary roddyExecutionStore1 roddyExecutionStore2 ...
#
# Compiles all configurations from all provided roddyExecutionStore directories.
# 
# * Output is JSON in which the simplified/essential configuration of the workflow for a set of directories
#   is displayed (=contexts).
# * If all .parameter files of an exec_* directory have the same configuration, only the exec_* directory
#   is reported in the "contexts" field.
# * If all exec_* directories in a roddyExecutionStore have the same configuration, only the
#   roddyExecutionStore is reported in the "contexts" field.
#
# Run the unit-tests with `pytest config-summary.py`
#
# Requirements: python2, more_itertools, pytest (for testing)
#
# Disclaimer: The plugin underwent quite some changes over the years. Use at your own responsibility.
#             If you find incorrect reportings, please file a bug report, such that the script gets
#             improved.
#
# TODO Have a script per workflow, to compress the parameters plus additional scripts in Roddy
#      to compile over many workflows.
# TODO Plugin mechanism for different workflows
# TODO Make plugin for other workflows
# TODO Diffing reports between groups to highlight differences for the user.
#
from __future__ import annotations

import sys
from collections import ChainMap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple, Union, Optional, TextIO, Mapping, Callable

import json
import logging
import re
from abc import ABCMeta, abstractmethod
from io import StringIO
from more_itertools import flatten

logger = logging.Logger(__name__)


def parse_parameter_line(line: str) -> Union[Tuple[str, str], None]:
    """
    Parse a line from Roddy's .parameter files.
    """
    match = re.match(r"declare(?:\s-\S)*\s+([^=]+)=(.+)", line.rstrip())
    if match:
        return match.groups()[0], match.groups()[1]
    else:
        return None


def test_parse_parameter_line():
    assert parse_parameter_line("if [[ -z \"${PS1-}\" ]]; then") is None
    assert parse_parameter_line("declare -x -i TARGETSIZE=74569526") == \
           ("TARGETSIZE", "74569526")
    assert parse_parameter_line("declare -x -i BWA_MEM_THREADS=8") == \
           ("BWA_MEM_THREADS", "8")
    assert parse_parameter_line("declare -x    SAMBAMBA_MARKDUP_VERSION=0.5.9") == \
           ("SAMBAMBA_MARKDUP_VERSION", "0.5.9")
    assert parse_parameter_line("declare -x    SAMBAMBA_MARKDUP_OPTS=\"-t 6 -l 9 " +
                                "--hash-table-size=2000000 --overflow-list-size=1000000 " +
                                "--io-buffer-size=64\"") == \
           ("SAMBAMBA_MARKDUP_OPTS",
            "\"-t 6 -l 9 --hash-table-size=2000000 --overflow-list-size=1000000 --io-buffer-size=64\"")
    assert parse_parameter_line("\n") is None


def parse_parameter_file(file: Path) -> Mapping[str, str]:
    """
    Parse the full parameter file and return a dictionary of variable names and values.
    Will return an empty dictionary, if no parameters can be parsed.
    """
    result = {}
    with open(file, "r") as f:
        for line in f.readlines():
            parsed = parse_parameter_line(line)
            if parsed is not None:
                result[parsed[0]] = parsed[1]
    return result


def read_version_info(input: TextIO) -> Mapping[str, str]:
    """
    Read a version info file and return a dictionary with the component (Roddy, plugins) and their
    versions parsed from the file.
    """
    result = {}
    first_line = input.readline().rstrip()
    match = re.match(r"Roddy version: (\d+\.\d+\.\d+|develop)", first_line)
    if not match:
        raise RuntimeError(f"Couldn't parse Roddy version from '{first_line}'")
    result["Roddy"] = match.groups()[0]
    input.readline()    # Drop "Library info:"
    for line in map(lambda l: l.rstrip(), input.readlines()):
        if line != "":
            match = re.match(r"Loaded plugin (\S+):(\S+)\sfrom.+", line)
            if not match:
                logger.error(f"Could not match plugin in: '{line}'")
            else:
                result[match.groups()[0]] = match.groups()[1]
    return result


def test_read_version_info():
    buffer = StringIO(
        """Roddy version: 3.5.9
Library info:
Loaded plugin AlignmentAndQCWorkflows:1.2.51-1 from ((/tbi/software/x86_64/otp/roddy/plugins/3.5/AlignmentAndQCWorkflows_1.2.51-1))
Loaded plugin COWorkflows:1.2.66-1 from ((/tbi/software/x86_64/otp/roddy/plugins/3.5/COWorkflows_1.2.66-1))
Loaded plugin PluginBase:1.2.1-0 from ((/tbi/software/x86_64/otp/roddy/roddy/3.5.9/dist/plugins/PluginBase_1.2.1))
Loaded plugin DefaultPlugin:1.2.2-0 from ((/tbi/software/x86_64/otp/roddy/roddy/3.5.9/dist/plugins/DefaultPlugin_1.2.2))
        """)
    parsed = read_version_info(buffer)
    assert parsed == {
        "Roddy": "3.5.9",
        "AlignmentAndQCWorkflows": "1.2.51-1",
        "COWorkflows": "1.2.66-1",
        "PluginBase": "1.2.1-0",
        "DefaultPlugin": "1.2.2-0"
    }


@dataclass
class ParameterSet:
    """
    The contexts are simply a list of directories for which the given parameters are valid.
    The client code may choose to simplify the context to a super-directory, to reduce the
    complexity of the output.
    """
    contexts: List[Path]
    parameters: Optional[Mapping[str, Mapping[str, Optional[str]]]]


class JobParameterEncoder(json.JSONEncoder):
    """
    An encoder to convert a JobParameter into a nested dict/list structure that can be serialized.
    """
    def default(self, obj):
        if isinstance(obj, ParameterSet):
            return {
                    "contexts": list(map(str, obj.contexts)),
                    "parameters": obj.parameters
            }
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


@dataclass
class ExecutionStore:
    """
    A representation of a roddyExecutionStore and the information logged in there.
    """
    roddy_store_path: Path
    execution_store_subdir: Path
    parameter_file_names: List[Path]
    parameters: List[ParameterSet] = field(default_factory=list)

    # The following methods return the path to the respective files/dirs starting with the
    # roddy_store_path.
    @property
    def execution_store(self) -> Path:
        return self.roddy_store_path / self.execution_store_subdir

    @property
    def versions_info_file(self) -> Path:
        return self.execution_store / "versionsInfo.txt"

    @property
    def parameter_files(self) -> List[Path]:
        return list(map(lambda p: self.execution_store / p,
                        self.parameter_file_names))

    def __repr__(self):
        result = [f"Execution store = {str(self.execution_store)}",
                  f"Version infos = {str(self.versions_info_file)}",
                  "Parameter files ="]
        for f in self.parameter_files:
            result.append(f"\t{str(f)}")

    # Factory methods and supplementary code.
    @classmethod
    def _collect_parameter_files(cls, execution_store: Path) -> List[Path]:
        return list(map(lambda f: f.relative_to(execution_store),
                        filter(lambda p: p.is_file() and p.name.endswith(".parameters"),
                               execution_store.iterdir())))

    @classmethod
    def from_path(cls, roddy_store_dir: Path, execution_store_subdir: Path) -> ExecutionStore:
        return ExecutionStore(
            roddy_store_path=roddy_store_dir,
            execution_store_subdir=execution_store_subdir,
            parameter_file_names=cls._collect_parameter_files(roddy_store_dir /
                                                              execution_store_subdir))


@dataclass
class RoddyStore:
    base_dir: Path
    execution_stores: List[ExecutionStore]
    # Whether there is a list or a single JobParameters object depends on whether the values are
    # consistent.
    parameters: List[ParameterSet] = field(default_factory=list)

    def __repr__(self):
        result = [f"Base directory = {self.base_dir}"]
        for store in self.execution_stores:
            result += store.__repr__()
        return "\n\n".join(result)

    # Factory methods and supplementary code.
    @classmethod
    def _collect_execution_dirnames(cls, roddy_store_path: Path) -> List[Path]:
        """
        Return dir names of the exec_ directories in the roddy_store_path.
        """
        return list(map(lambda p: Path(p.name),
                        filter(lambda p: ((p.is_dir() or p.is_symlink()) and
                                          p.name.startswith("exec_")),
                               roddy_store_path.iterdir())))

    @classmethod
    def from_path(cls, base_dir: Path) -> RoddyStore:
        stores = []
        for exec_store_subdir in cls._collect_execution_dirnames(base_dir):
            stores.append(ExecutionStore.from_path(base_dir, exec_store_subdir))
        return RoddyStore(base_dir=base_dir,
                          execution_stores=stores)


def subdict_str(the_dict: Mapping[str, str], keys: List[str], default=None) \
        -> Mapping[str, Optional[str]]:
    """
    Just a helper function that allows subsetting the keys of a Dict[str,str].
    Note that key that are requested but not in the input dictionary will be returned with
    value Nane.
    """
    return {k: the_dict.get(k, default) for k in keys}


def simple_combine_parameter_sets(parameters: List[ParameterSet]) \
        -> List[ParameterSet]:
    """
    This is a strategy for combining parameter sets. It's very simple: Just check whether all
    parameters are identical and return one representative, if they are. Otherwise, return the
    original input.
    """
    first_parameters = parameters[0].parameters
    consistent = all(map(lambda job: job.parameters == first_parameters,
                         parameters))
    if consistent:
        return [ParameterSet(contexts=list(flatten(map(lambda p: p.contexts,
                                                       parameters))),
                             parameters=first_parameters)]
    else:
        return parameters


def test_simple_combine_parameter_sets():
    # Two equal
    assert simple_combine_parameter_sets([ParameterSet(contexts=[Path("b")],
                                                       parameters={"b": 1}),
                                          ParameterSet(contexts=[Path("c")],
                                                       parameters={"b": 1})]) == \
        [ParameterSet(contexts=[Path("b"), Path("c")],
                      parameters={"b": 1})]

    # Two unequal
    assert simple_combine_parameter_sets([ParameterSet(contexts=[Path("b")],
                                                       parameters={"b": 2}),
                                          ParameterSet(contexts=[Path("c")],
                                                       parameters={"b": 1})]) == \
           [ParameterSet(contexts=[Path("b")],
                         parameters={"b": 2}),
            ParameterSet(contexts=[Path("c")],
                         parameters={"b": 1})]

    # Three, two equal, one not
    assert simple_combine_parameter_sets([ParameterSet(contexts=[Path("b")],
                                                       parameters={"b": 2}),
                                          ParameterSet(contexts=[Path("c")],
                                                       parameters={"b": 1}),
                                          ParameterSet(contexts=[Path("d")],
                                                       parameters={"b": 2})]) == \
           [ParameterSet(contexts=[Path("b")],
                         parameters={"b": 2}),
            ParameterSet(contexts=[Path("c")],
                         parameters={"b": 1}),
            ParameterSet(contexts=[Path("d")],     # This one is not grouped with the first!
                         parameters={"b": 2})]


def grouping_combine_parameter_sets(job_parameters: List[ParameterSet]) \
        -> List[ParameterSet]:
    """
    This strategy for combining results, makes an all vs. all comparison of parameters and groups
    those that are identical. The contexts (usually exec_* dirs) are compiled in a list.
    """
    # We use the string-mapped parameters, because we want to parameters that are *identical*.
    mapped_params = map(lambda jp: (json.dumps(jp.parameters, sort_keys=True), jp),
                        job_parameters)

    params_by_params: dict = {}
    for key, params in mapped_params:
        if key in params_by_params.keys():
            # Accumulate the contexts for these parameters.
            params_by_params[key] = {
                "contexts": params_by_params[key]["contexts"] + params.contexts,
                "parameters": params_by_params[key]["parameters"]
            }
        else:
            # Initialize the context for these parameters.
            params_by_params[key] = {
                "contexts": params.contexts,
                "parameters": params.parameters
            }

    result = list(map(lambda p: ParameterSet(**p), params_by_params.values()))
    return result


def test_grouping_combine_parameter_sets():
    # Two equal
    assert grouping_combine_parameter_sets([ParameterSet(contexts=[Path("b")],
                                                         parameters=dict({"b": 1})),
                                            ParameterSet(contexts=[Path("c")],
                                                         parameters=dict({"b": 1}))]) == \
           [ParameterSet(contexts=[Path("b"), Path("c")],
                         parameters=dict({"b": 1}))]

    # Two unequal
    assert grouping_combine_parameter_sets([ParameterSet(contexts=[Path("b")],
                                                         parameters={"b": 2}),
                                            ParameterSet(contexts=[Path("c")],
                                                         parameters={"b": 1})]) == \
           [ParameterSet(contexts=[Path("b")],
                         parameters={"b": 2}),
            ParameterSet(contexts=[Path("c")],
                         parameters={"b": 1})]

    # Three, two equal, one not
    assert grouping_combine_parameter_sets([ParameterSet(contexts=[Path("b")],
                                                         parameters={"b": 2}),
                                            ParameterSet(contexts=[Path("c")],
                                                         parameters={"b": 1}),
                                            ParameterSet(contexts=[Path("d")],     # same paras as b!
                                                         parameters={"b": 2})]) == \
           [ParameterSet(contexts=[Path("b"), Path("d")],    # b and d are grouped
                         parameters={"b": 2}),
            ParameterSet(contexts=[Path("c")],
                         parameters={"b": 1})]


class VersionAnalysisStrategy(metaclass=ABCMeta):

    @abstractmethod
    def run(self, roddy_store: RoddyStore) -> Union[ParameterSet, List[ParameterSet]]:
        """
        Multiple execution stores may exist for an analysis directory, if the workflow was run multiple
        times. This function tries to combine the information from multiple such runs. If the
        information cannot be reconciled. The rules are as follows

        * Versions in later runs override versions in earlier runs.
        * Cluster jobs are generally reported separately, unless they all used the same versions.
        """
        pass


class AlignmentVersionAnalysisStrategy(VersionAnalysisStrategy):
    """
    This encodes the interpretation of the alignment configuration and mapping to parameter sets.
    Furthermore, (and this could thus be extracted into its own class) it combines the results
    using a combination strategy.
    """

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
        if parameters.get(self.runFastQC, "false") == "true" \
                and parameters.get(self.runAlignmentOnly, "false") == "false":
            return subdict_str(parameters,
                               [self.runFastQC,
                                self.runFastQCOnly,
                                self.FASTQC,
                                self.JAVA])
        else:
            return {}

    def trimmomatic_parameters(self,
                               plugin_versions: Mapping[str, str],
                               parameters: Mapping[str, str]) -> Mapping[str, Optional[str]]:
        if parameters.get(self.useAdapterTrimming, "false") == "true":
            return {
                "TRIMMOMATIC_VERSION": "0.30",  # Hardcoded!
                **subdict_str(parameters,
                              [self.useAdapterTrimming,
                               self.ADAPTOR_TRIMMING_OPTIONS_0,
                               self.ADAPTOR_TRIMMING_OPTIONS_1,
                               self.CLIP_INDEX,
                               self.JAVA])
            }
        else:
            return {}

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

    def _interpret_parameters(self,
                              plugin_versions: Mapping[str, str],
                              parameters: Mapping[str, str]) -> \
            Mapping[str, Mapping[str, Optional[str]]]:
        """
        This function recapitulates the logic of the workflow. Only parameters and versions are
        returned that were actually relevant for the analysis described by the parameters and
        plugin_versions. Therefore, version-specific differences may also be addressed here.

        Returns a nested dictionary structure with string values, that describe the parameters.
        """
        # Some parameters may not be set and would otherwise be reported as "null". Their defaults
        # are generally "false". We fill these values in, in case they are missing.
        parameters = ChainMap(parameters, self.possibly_missing_value_defaults)
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
            "wes": self.wes_parameters(plugin_versions, parameters),
            "wgbs": self.wgbs_parameters(plugin_versions, parameters),
        }
        return results

    def run(self, roddy_store: RoddyStore) -> List[ParameterSet]:
        """
        Collect version and other run information from all execution store directories. For each
        exec_ directory, return a dictionary of key value pairs with information. Information is
        left out, if they are not relevant. E.g. a software version may be omitted, if the software
        was not used, because it was turned off some other parameter value.
        """
        per_exec_store: List[ParameterSet] = []
        for exec_store in roddy_store.execution_stores:
            with open(exec_store.versions_info_file, "r") as f:
                versions = read_version_info(f)

            if len(exec_store.parameter_files) == 0:
                store_parameters = [ParameterSet(contexts=[roddy_store.base_dir],
                                                 parameters={})]
            else:
                job_parameters: List[ParameterSet] = []
                for job_parameters_file in exec_store.parameter_files:
                    raw_parameters = parse_parameter_file(job_parameters_file)
                    job_parameters.append(ParameterSet([job_parameters_file],
                                                       self._interpret_parameters(versions,
                                                                                  raw_parameters)))

                # In principle, each job may have different parameters. We first combine the parameters
                # of all .parameter files. This is basically done by rejecting the exec_* directory
                # because finding divergent parameters is a serious sign of messing with Roddy, which
                # we don't support.
                store_parameters = grouping_combine_parameter_sets(job_parameters)
                if len(store_parameters) > 1:
                    # We issue an error, if different parameters are used for different workflow runs
                    # on the same output directory, because this is a serious indication of messing
                    # with Roddy's execution. But we do not stop the whole analysis.
                    print(f"Inconsistent job parameters in {exec_store.execution_store}.",
                          file=sys.stderr)
                else:
                    # To simplify the output, we replace the list of parameter files by the execution
                    # store directory. This expresses also to the client, that all execution store
                    # directories contain the same parameter set.
                    store_parameters = [ParameterSet(contexts=[Path(exec_store.execution_store)],
                                                     parameters=store_parameters[0].parameters)]

            exec_store.parameters = store_parameters
            per_exec_store += store_parameters

        # What is much more likely is that the configurations from different executions of Roddy,
        # differ, even though, usually, all jobs will have the same parameters. We try to combine
        # all parameters from all executions, but if that is not possible, report all variants.
        roddy_store_parameters = grouping_combine_parameter_sets(per_exec_store)

        # Like before, we look at all execution stores, and if the all have the same parameters,
        # we can simplify to only list the roddyExecutionStore directory as context.
        if len(roddy_store_parameters) == 1:
            roddy_store_parameters = [ParameterSet(contexts=[roddy_store.base_dir],
                                                   parameters=roddy_store_parameters[0].parameters)]

        roddy_store.parameters = roddy_store_parameters
        return roddy_store_parameters


if __name__ == "__main__":
    """
    Just call the script with a set of `roddyExecutionStore` directories, collect the version 
    information and generate a report.
    
    If all directories and recursively, exec_* directories and .paramater files used the same
    configuration, only this config is returned in an interpreted way that minimizes the output.
    
    If there are exec_* directories with different configurations in the same roddyExecutionStore
    then they will be reported separately (the configurations are inconsistent).
    
    If there are even .parameter files in the same exec_* directory, you messed up a single run
    of Roddy. S.b. probably messed around with the configs. Pray that the results are usable!
    
    DISCLAIMER: This script can only find version information as used by the workflow itself. If 
                you hacked and ran partial analyses manually with other tools, then this script
                will obviously produce wrong results.
    """
    directories = sys.argv[1:]
    if len(directories) == 0:
        print("No roddyExecutionStore directories provided", file=sys.stderr)
        sys.exit(1)
    analyser = AlignmentVersionAnalysisStrategy()
    per_roddy_store = []
    for roddy_execution_store_dir in directories:
        exec_store = RoddyStore.from_path(Path(roddy_execution_store_dir))
        per_roddy_store += analyser.run(exec_store)

    combined: List[ParameterSet] = grouping_combine_parameter_sets(per_roddy_store)
    if len(combined) > 1:
        print("Could not combine configurations in roddyExecutionStores", file=sys.stderr)
        json.dump(combined, sys.stdout, cls=JobParameterEncoder)
        sys.exit(2)
    else:
        print("Same configuration in roddyExecutionStores", file=sys.stderr)
        json.dump(combined, sys.stdout, cls=JobParameterEncoder)
        sys.exit(0)


