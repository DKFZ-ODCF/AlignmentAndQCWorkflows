# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/TheRoddyWMS/Roddy) workflows.

This is the branch for 1.2.73 maintenance branch with only minimal documentation. Please see the master branch for details.

## Read-Reordering with Unidirectional WGBS Data

By setting `reorderUnidirectionalWGBSReadPairs` the read-reordering script will be run that decides based on the relative frequencies of TC and AG dinucleotides in both reads, what is the most likely correct orientations of the reads, and may swap the two reads.

Note that after the swapping, the read-numbers are reversed. What was R1 in the input FASTQ will be R2 in the output BAM, and vice versa.

Furthermore, not all reads can be unambiguously classified. These unclassified reads are currently dropped. 

The original script with a documentation of the underlying ideas can be found [here](https://github.com/cimbusch/TWGBS.git).

## Change Logs

* 1.2.73-3 (branch-specific change)
  - Updated unidirectional WGBS read-reordering script from [here](https://github.com/cimbusch/TWGBS.git)
  - Actually include the WGBS read-reordering script

* 1.2.73-2 (branch-specific changes)
  - Improved error checking and reporting for BWA and surrounding pipe
  
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
  - Fingerprinting (not for WGBS, yet)
  - Bugfix affecting CLIP_INDEX in configuration 
