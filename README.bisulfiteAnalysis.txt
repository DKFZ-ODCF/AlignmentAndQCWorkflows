== Description

Whole Genome Bisulfite Sequencing (WGBS) workflow based on methylCtool.

== Run flags / switches

Switch                      Default Description


== Changelist

* Version update 1.2.73

  - lifted to Roddy 2.4
  - support for loading environment modules (via DefaultPlugin 1.2.0)
  - removed references to PBS to allow running on LSF cluster
  - fixed FASTQC code
  - fingerprinting


* Version update 1.1.39
 
  - Initial WGBS support. Duplication marking with Picard or Sambamba.
