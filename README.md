# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/eilslabs/Roddy) workflows:

- PanCancer alignment workflow for whole genome (WGS) and exome (WES)
- Bisulfite core workflow (WGBS) using [methylCtools](https://github.com/hovestadt/methylCtools)

These are basically BWA alignment (bwa mem) workflows with plenty of additional quality control steps. QC-steps acting on the BAM files are mostly done through piping the input data into multiple QC tools simultaneously to achieve a high performance. All workflows can be run for each sample or with combinations of tumor- and control-samples.

## Documentation

Please refer to the [wiki](../../wiki/) for detailed information on all aspects of the workflow, including installation, configuration structure and interpretation of results.  

## German Network on Bioinformatic Infrastructure (de.NBI) 

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=AlignmentAndQCWorkflows">this very short (30-60 seconds) survey</a> you support our efforts to improve this tool.</td></tr></table>

## Change Logs

Please see the [changelogs file](Changelog.md) for details.
