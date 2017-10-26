== Description

The workflow is based on the original qc genome workflow and uses most of its run flags.
It utilizes a set of optimized methods leading to a very reduced workflow with only three steps.
In addition it allows the usage of externally set bam files and fastq files.

== Run flags / switches

Switch                      Default Description
runCoveragePlotsOnly        false

Value                       Description
overrideFastqFiles          A list of
overrideSamples
overrideBamFiles


== Changelist

* Version update 1.2.51

  - lifted to Roddy 2.4
  - support for loading environment modules (via DefaultPlugin 1.2.0)
  - removed references to PBS to allow running on LSF cluster
  - fixed FASTQC code

...

- Renamed plugin to AlignmentAndQCWorkflows.

* Version update to 1.0.186

- Sambamba support for duplication marking.

* Version update to 1.0.182

- The resource requests for jobs were increased to be on the safe side when switching to the new cluster (in particular concerning OOM killer).
- Adapted QCWF scripts to match PBS_QUEUE on convey* rather than "convey", to ensure operation on the new cluster, where queues are named convey_fast, convey_medium and convey_long.

* Version update to 1.0.180

* Version update to 1.0.178

* Version update to 1.0.177

- Calculate MD5 sums for both, Picard and biobambam-based
  workflows, using md5sum in a separate branch of pipes.
- Error-checks after mv commands in alignAndPairSlim.sh
  and mergeAndMarkOrRemoveDuplicatesSlim.sh.

* Version update to 1.0.173

- qualitycontrol.json

* Version update to 1.0.168

- The biobambam branch of the slim mark duplicates script
  (mergeAndMarkOrRemoveDuplicatesSlim.sh) now produces
  merged BAM md5sum file.

* Version update to 1.0.166

- Removed requests for the lsdf from all scripts.

* Version update to 1.0.164

- Added fastq_list configuration value that allows to override inputDirectories
  and directly provide FASTQs on the commandline via --cvalues.

* Version update to 1.0.161

* Version update to 1.0.158

* Version update to 1.0.135

* Version update to 1.0.132

* Version update to 1.0.131

* Version update to 1.0.114

* Version update to 1.0.109

* Version update to 1.0.105

* Version update to 1.0.104

* Version update to 1.0.103
