# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/TheRoddyWMS/Roddy) workflows.

This is the branch for 1.2.73 maintenance branch with only minimal documentation. Please see the master branch for details.

## Read-Reordering with Unidirectional WGBS Data

By setting `reorderUnidirectionalWGBSReadPairs` the read-reordering script will be run that decides based on the relative frequencies of TC and AG dinucleotides in both reads, what is the most likely correct orientations of the reads, and may swap the two reads.

Note that after the swapping, the read-numbers are reversed. What was R1 in the input FASTQ will be R2 in the output BAM, and vice versa.

Furthermore, not all reads can be unambiguously classified. These unclassified reads are currently dropped. 

The original script with a documentation of the underlying ideas can be found [here](https://github.com/cimbusch/TWGBS.git).

## Change Logs

* 1.2.73-203 (branch-specific change)
  - minor: Optional ALT-chromosome processing via bwa.kit's `bwa-postaln.js`.
    * `ALT_FILE`: Defaults to be `$INDEX_PREFIX.aln`
    * `K8_VERSION`: Used by `tbi-lsf-cluster.sh` environment script.
    * `K8_BINARY`: Needs to be set, e.g. in the environment setup script. Done for `tbi-lsf-cluster.sh` based on `K8_VERSION`.
    * Set `runBwaPostAltJs=true` to activate the ALT chromosome processing. Default: `false`.
    * Set `bwaPostAltJsPath` to point to the `bwa-postalt.js` script.
    * You can set `bwaPostAltJsParameters`.
  - minor: Set `useCombinedAlignAndSampe=false` and `runSlimWorkflow=true` in the default config. The other option is unmaintained and obviously wasn't used for years.
  - minor: Set `workflowEnvironmentScript` to `workflowEnvironment_tbiLsf`. The previous value was reasonable only as long our PBS cluster existed. The `tbi-pbs-cluster.sh` script is also removed.
  - patch: Set some resources limits for "fastqc" job.

* 1.2.73-202 (branch-specific change)
  - patch: call of `grDevices.pdf()` in `chrom_diff.r` lead to unexpected 
    (truncated to at most 511 characters) file name of output pdf 
  
* 1.2.73-201 (branch-specific change)
  - minor: explicit `MBUFFER_VERSION` can be specified for tbi-lsf environment
  
* 1.2.73-2, 1.2.73-200 (branch-specific change)
  - minor: Included unidirectional WGBS read-reordering script from [here](https://github.com/cimbusch/TWGBS.git)
  - patch: Improved error checking and reporting for BWA and surrounding pipe
  
* 1.2.73-1 (branch-specific changes)
  - major: Lifted to Roddy 3.0 release (official LSF-capable release)
  - patch: Bugfix with wrong Bash function export

* 1.2.73
  - major: Lifted 1.1.73 to Roddy 2.4 (development-only release)
  - minor: Fingerprinting support also for WGBS
  - minor: sambamba 0.5.9 for sorting and viewing BAMS
  - patch: BAM termination sequence check

* 1.1.73
  - minor: Fingerprinting (not for WGBS, yet)
  - patch: Bugfix mergeOnly step WGBS
  - patch: Substituted sambamba-based compression by samtools compression for improved stability, time, and memory consumption
  - patch: Tuning (tee -> mbuffer)
  - patch: Node-local scratch by default
  - patch: Bugfix affecting CLIP_INDEX in configuration 
