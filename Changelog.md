# Versioning Policy

## Release Branches

Various versions of the workflow are or have been in production mode at the DKFZ/ODCF as part of the data managment system [OTP](https://otp.dkfz.de/otp/). The production versions often have dedicated release branches named "ReleaseBranch_$major.$minor\[.$patch\]" in which a only restricted changes have been made in order to ensure compatibility of the results:

  * No changes that alter previous output.
  * New important features are sometimes backported -- as long as they do not change the previous results.
  * Bugfixes that allow running the workflow on some data on which it previously crashed, but that do not alter the existing output, are included.

## Roddy 2 to 3 Shift

For the shift from Roddy 2 to Roddy 3 workflow version were increased in the minor number. Specifically this happened to the following versions:

  - 1.1.51 -> 1.2.51
  - 1.1.73 -> 1.2.73
  - 1.0.182 -> 1.2.182
  
> Therefore note that [ReleaseBranch_1.2.182](../../tree/ReleaseBranch_1.0.182) is __not the newest branch, but the oldest__! It was derived from a very old version of the workflow ([QualityControlWorkflows_1.0.182](../../tree/ReleaseBranch_1.0.182)) at a time where the versioning system was not fixed to [semver 2.0](https://semver.org/).

## Semantic Versioning 2.0

Starting with version 1.2.76 we switched to [Semantic Versioning 2.0](https://semver.org/) with a focus on user-oriented changes. This means the version numbers are increased according to semantic **changes on the interface to the user**, that is variables and output. The compatibility management with Roddy and upstream plugin versions is automatically managed. 

In exception to this strategy backports etc. for maintenance branches are created by suffixing a number separated by '-' to the semantic version.

# Change Logs

* 1.4.0
  - Recognize truncated FASTQs from BWA error log
  - Update to Roddy 3.5
  - Major documentation update (workflow structure plots)
  - JVM code refactorings
  - Merged in changes from 1.2.73-2
    - Updated undirectional read-reordering script and integrated into WGBS pipeline
    - Improved error checking and reporting for BWA and surrounding pipe (set -e, PID registry)
    - Got (most of) the BWA methyl-seq code to run with set -e to improve error robustness and handling.
    - Imported BamToFastqPlugin tempfile and process ID registry code (already well tested)
    - Backported bash unit tests from master
    - Bash pipe extension framework

* 1.3.0
  - Coverage separate for mouse and human in xenograft assemblies
  - Check ACEseq QC input file before start
  - Use externally provided trimmomatic
  - Update to Roddy 3.2 dependency
  - Better documentation
  - Conda environment
  - Further removal of unused/unmaintained scripts
  - Better tolerance to small datasets (mostly for testing)
  - Fixed some Perl tests
  - Diverse bugfixes
  - MIT licence, where necessary copyleft licence

* 1.2.76
  - Lifted 1.1.76 to Roddy 3
  - LSF support
  - Support for loading environment modules (via DefaultPlugin 1.2.2)
  - Check input BAMs for syntactic completeness (BAM trailer)
  - Check FASTQ files before submission
  - Classify FASTQs as QC-passed or QC-failed based on FASTQC output
  - ACEseq QC (runACEseqQc:Boolean, GC_CONTENT_FILE_ALN, REPLICATION_TIME_FILE_ALN, MAPPABILITY_FILE_ALN, CHROMOSOME_LENGTH_FILE_ALN); not on WES
  - Stabilization WGBS
  - Fingerprinting for WGBS
  - Turned off faulty fingerprinting on Conveys
  - Additionally to upper-case SAMPLE, RUN, etc. also send lower-case versions to jobs  
  - Refactorings and code cleanup
  - Deleted old BWA sampe code
  - Renamed alreadyMergedLanes.pl to missingReadGroups.pl
  - Removed bwaErrorCheckingScript; now function in workflowLib.sh
  - Documentation (Readme.md)
  - Unit tests for Bash functions in workflowLib.sh and bashLib.sh (based on [shunit2](https://github.com/kward/shunit2.git))
  - Refactoring: bamFileExists -> useOnlyExistingTargetBam
  - bugfix: FASTQC code
  - bugfix: use existing BAM files
  - bugfix: BWA error recognition

* 1.2.73-2 (branch-specific changes)
  - Updated undirectional read-reordering script and integrated into WGBS pipeline
  - Improved error checking and reporting for BWA and surrounding pipe (set -e; PID registry)
  - Got (most of) the BWA methyl-seq code to run with set -e to improve error robustness and handling.
  - Imported BamToFastqPlugin tempfile and process ID registry code (already well tested)
  - Backported bash unit tests from master
  - Bash pipe extension framework
   
* 1.2.73-1 (branch-specific changes)
  - Lifted to Roddy 3.0 release (official LSF-capable release)
  - Bugfix with wrong Bash function export

* 1.2.73
  - Lifted 1.1.73 to Roddy 2.4 (development-only release)
  - Fingerprinting support also for WGBS
  - sambamba 0.5.9 for sorting and viewing BAMS
  - BAM termination sequence check

* 1.1.73
  - Bugfix mergeOnly step WGBS
  - Substituted sambamba-based compression by samtools compression for improved stability, time, and memory consumption
  - Tuning (tee -> mbuffer)
  - Node-local scratch by default
  - Fingerprinting for WES and WGS (runFingerprinting:Boolean, fingerprintingSitesFile); not for WGBS yet
  - Bugfix affecting CLIP_INDEX in configuration 
  - Tuned parameters for sambamba support and extracted BAM compression into separate step for performance reasons

* 1.2.51-2 (branch-specific changes)
  - Improved error checking and reporting for BWA and surrounding pipe

* 1.2.51-1 (branch-specific changes)
  - Update to Roddy 3.0 release (official LSF-capable release)
  - Bugfix in tbi-lsf-cluster.sh

* 1.2.51
  - Lifted 1.1.51 to Roddy 2.4 (development-only release)
  - FASTQ quality classification (Xavier)
  - BAM termination sequence check
  - Bugfixes in WGBS (off-by-one, meth-call splitting CG/CH)
  - Further bugfixes

* 1.1.51
  - Improved and error checking
  - Pre-submission executability checks
  - Tuning (sambamba flagstat -t 1)
  - Use local scratch on nodes
  - Resource size 't' for testing purposes
  - Progress on WGBS workflow
  - Fixed off-by-one error in moabs output

* Version upgrade to 1.1.39  
 
- Initial WGBS support. Duplication marking with Picard or Sambamba.

* 1.1.2
  - Rename from QualityControlWorkflows to AlignmentAndQCWorkflows
  - First development version of WGBS workflow
  - Adaptation to Roddy API changes (FileSystemAccessProvider)


## QualityControlWorkflows

* 1.0.186
  - Updated dependency on COWorkflows version 1.1.23
  - Sambamba support for duplication marking
  - Fixed the merging of new lanes into existing merged BAM
  - Generalized the flagstats parser (samtools 1+ format with "supplementary reads")
  - Refactorings
  - workflowLib.sh for workflow-specific Bash code
  - Increased resource requirements  
  - Resource size 't' (for testing purposes)
  - Flagstats support for samtools 1.0+ and sambamba 0.5.9 (supplementary reads)

* 1.2.182 (branch-specific changes)
  - Roddy 3 support (official LSF-capable release)
  - BAM termination sequence check

* 1.0.182-1
  - Increased resource requirements
  - `bam` configuration value to provide an externally located BAM file as initial merged BAM into which to merge additional lane-BAMs
  - Fixed missing metric file required for QC summary file creation issue
  - Bugfixes

* 1.0.182
  - Increased resource requirements
  - Refactoring
  - Adapted QCWF scripts to match PBS_QUEUE on convey* rather than "convey" (to match "convey*" queue names)

* 1.0.180
  - Imported compiled coverageQc binary
  - Added tiny/testing (t) resource set for WGS alignment workflow
  - Added qcJson.pl call to targetExtractCoverageSlim job
  - Per-Read Group post merge QC added
  - Git repo created from original SVN checkout
  - WES: Added qualitycontrol_targetExtract.json

* 1.0.178

* 1.0.177
  - Calculate MD5 sums for both, Picard and Biobambam-based workflows, using md5sum in a separate branch of pipes.
  - Error-checks after `mv` commands in alignAndPairSlim.sh and mergeAndMarkOrRemoveDuplicatesSlim.sh.

* 1.0.173
  - Added qualitycontrol.json

* 1.0.168
  - The biobambam branch of the slim mark duplicates script (mergeAndMarkOrRemoveDuplicatesSlim.sh) now produces merged BAM md5sum file.

* 1.0.166
  - Removed requests for the lsdf from all scripts.

* 1.0.164
  - Added fastq_list configuration value that allows to override inputDirectories and directly provide FASTQs on the commandline via --cvalues.

* 1.0.161

* 1.0.158

* 1.0.135

* 1.0.132

* 1.0.131

* 1.0.114

* 1.0.109

* 1.0.105

* 1.0.104

* 1.0.103
