== Description

Please take a look at the genome quality control readme file.
The workflows are the same.

== Changelist

* Version update to 1.2.76

- lifted to Roddy 3.0
- support for loading environment modules (via DefaultPlugin 1.2.0)
- removed references to PBS to allow running on LSF cluster
- fixed FASTQC code
- ACEseq QC (runACEseqQc:Boolean, GC_CONTENT_FILE_ALN, REPLICATION_TIME_FILE_ALN, MAPPABILITY_FILE_ALN, CHROMOSOME_LENGTH_FILE_ALN)
- fingerprinting for WGBS
- turned off faulty fingerprinting on Conveys
- classify FASTQs as QC-passed or QC-failed based on FASTQC output

* Version update to 1.1.73

- fingerprinting for WES and WGS (runFingerprinting:Boolean, fingerprintingSitesFile)
- tuned parameters for sambamba support and extracted BAM compression into separate step for performance reasons

* Version update to 1.1.51

- sambamba support
- use local scratch on nodes ()
- resource size 't' for testing purposes
- fixed off-by-one error in moabs output

* AlignmentAndQCWorkflows-1.1.2

- Renamed plugin to AlignmentAndQCWorkflows.

* Version update to 1.0.186

- sambamba support for duplication marking
- resource size 't' for testing purposes
- support for samtools 1.0+ and sambamba 0.5.9 (supplementary reads)

* Version update to 1.0.182

* Version update to 1.0.180

- qualitycontrol_targetExtract.json

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
