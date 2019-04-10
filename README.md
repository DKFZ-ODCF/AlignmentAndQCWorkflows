# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/eilslabs/Roddy) workflows:

- PanCancer alignment workflow for whole genome (WGS) and exome (WES)
- Bisulfite core workflow (WGBS) using [methylCtools](https://github.com/hovestadt/methylCtools)

These are basically BWA alignment (bwa mem) workflows with plenty of additional quality control steps. QC-steps acting on the BAM files are mostly done through piping the input data into multiple QC tools simultaneously to achieve a high performance. All workflows can be run for each sample or with combinations of tumor- and control-samples.

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=AlignmentAndQCWorkflows">this very short (30-60 seconds) survey</a> you support our efforts to improve this tool.</td></tr></table>

Please refer to the [wiki](wiki/) for detailed information on all aspects of the workflow, including installation, configuration structure and interpretation of results.  

# Release Branches

Various versions are or have been in production mode at the DKFZ/ODCF. These often have dedicated release branches named "ReleaseBranch_$major.$minor\[.$patch\]" in which only certain changes have been made:

  * No changes that alter previous output.
  * New important features are sometimes backported -- as long as they do not change the previous results.
  * Bugfixes that allow running the workflow on some data on which it previously crashed, but that do not alter the existing output, are included.
  
> Note that [ReleaseBranch_1.2.182](../../tree/ReleaseBranch_1.0.182) is __not the newest branch, but the oldest__! It was derived from a very old version of the workflow ([QualityControlWorkflows_1.0.182](../../tree/ReleaseBranch_1.0.182)) at a time where the versioning system was not fixed to [semver 2.0](https://semver.org/).

# Change Logs

Please see the [changelogs file](Changelog.md) for details.