== Description

Some description

== Run flags / switches

Switch                      Default Description
runCoveragePlots            true    Create coverage plots at the end of the workflow
runFastQ                    true    Run the fastqc tool on fastq files

runAlignmentOnly            false   Run only the alignment steps
runCollectBamFileMetrics    false   Collect bam metrics for merged bams
useCombinedAlignAndSampe    false   Use bwa mem instead of bwa
useExistingPairedBams       false   Do not create paired bams and use existing ones. Merge will still be done!
runExomeAnalysis            false   Does the workflow run for genome or exome analysis
runFastQCOnly               false   Stop after the fastq tool
runSlimWorkflow             false   Use a full workflow with a lot of separate steps or a slim workflow with a lot of in-job piping.

== Changelist

* Version update to 1.2.51

* Version update to 1.0.51

- Renamed plugin to AlignmentAndQCWorkflows.

* Version update to 1.0.186

* Version update to 1.0.182

* Version update to 1.0.180

* Version update to 1.0.178

* Version update to 1.0.177

* Version update to 1.0.175

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
