# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/eilslabs/Roddy) workflows:

- PanCancer alignment workflow for whole genome (WGS) and exome (WES)
- Bisulfite core workflow (WGBS) using [methylCtools](https://github.com/hovestadt/methylCtools)

These are basically BWA alignment (bwa mem) workflows with plenty of additional quality control steps. QC-steps acting on the BAM files are mostly done through piping the input data into multiple QC tools simultaneously to achieve a high performance. All workflows can be run for each sample or with combinations of tumor- and control-samples.

# Configuration

All configuration variables are documented in the XML files in the `resources/configurationFiles/` directory. There is one XML for each workflow. Note that the workflows depend on each other, i.e. the WES and WGBS workflows are extension of the WGS workflow -- this can be recognized from the `imports` attribute of the top-level configuration tag in the XMLs. This means that most options of the WGS workflow also affect the other two workflows. Conversely, settings in the WGBS and WES workflow may override those in the WGS workflow. Some processing steps, notably those of the ACEseq quality control (QC) are not valid for the WES workflow. Note that the plugin depends on the [COWorkflowBasePlugin](https://github.com/TheRoddyWMS/COWorkflowsBasePlugin), which has its own configurations affecting this alignment plugin.

See the [Roddy](https://github.com/TheRoddyWMS/Roddy) documentation for a description of how to configure and run workflows.

# Software Requirements

## Conda

The workflow contains a description of a [Conda](https://conda.io/docs/) environment. A number of Conda packages from [BioConda](https://bioconda.github.io/index.html) are required. You should set up the Conda environment at a centralized position available from all compute hosts. 

First install the BioConda channels:
```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Then install the environment

```
conda env create -n AlignmentAndQCWorkflows -f $PATH_TO_PLUGIN_DIRECTORY/resources/analysisTools/qcAnalysis/environments/conda.yml
```

The name of the Conda environment is arbitrary but needs to be consistent with the `condaEnvironmentName` variable. The default for that variable is set in `resources/configurationFiles/qcAnalysis.xml`.

## Disclaimer for the Conda Version

The software versions in the Conda environment are yet not exactly the same as the ones in our local HPC infrastructure at the German Cancer Research Center. Among the many differences the following are probably most interesting for the user of the data:

|Package   | DKFZ version | Conda version | Comment                |
|----------|--------------|---------------|------------------------|
|biobambam | 0.0.148      | 2.0.79        | As long as you do not select `markDuplicatesVariant=biobambam` this won't be a problem, as biobambam is only used for sorting BAMs. Note further, we did not manage to get bamsort 2 from Conda to run on a CentOS 7 VM. You can also use `useBioBamBamSort=false` to sort with samtools.|
|picard    | 1.125        | 1.126         | Probably no big deal. |
|bwa       | patched 0.7.8| 0.7.8         | For the WGBS workflow we currently use a patched version of BWA that does not check for the "/1" and "/2" first and second read marks. This version is not available in BioConda and thus the WGBS workflow won't work with the Conda environment. |
|R         | 3.4.0        | 3.4.1         | Probably no big deal. |
|trimmomatic| 0.30        | 0.33          |                       | 

We successfully tested the Conda environment imported as described above and using the parameters `useBioBamBamSort=false`, `markDuplicatesVariant=sambamba`, `workflowEnvironmentScript=workflowEnvironment_conda` and `condaEnvironmentName=AlignmentAndQCWorkflows` on WGS data.

## WGBS Data Processing and methylCtools

The current implementation of the WGBS workflow uses [methylCtools](https://github.com/hovestadt/methylCtools) requires a patched BWA version. Note that the methylCtools version in this repository is not as new as the one in the official Github repository. Furthermore, 

## Recompiling the D-based Components

Two programs in this repository -- `genomeCoverage.d` and `coverageQc.d` -- were written in the programming language [D](https://dlang.org/) and are provided as binaries and source code. If the need arises to recompile them you can find the build instructions in `resources/`. For the compilation you will need

  * the D-compiler [LDC 0.12.1](https://github.com/ldc-developers/ldc/releases/tag/v0.12.1) compiler
  * and [BioD](https://github.com/lomereiter/BioD) master branch (commit 8b633de) 

# Resource Requirements

The workflow is rather tuned to minimize IO. For instance, the tools are glued together using pipes. However, the duplication marking and the BAM sorting steps produce temporary files. These two and the BWA step are also the memory-hungry steps, while BWA is the step that requires most CPU time. 
 
Have a look at the resource definitions in XMLs in the `resources/configurationFiles/` directory, which are rather conservative and will cover almost all 30x-80x human data sets we received from our X10 sequencers. The XMLs contain multiple parameters that allow you can tweak the actually used memory and cores. 

* bwa (BWA_MEM_THREADS)
* mbuffer (MBUFFER_SIZE_LARGE, MBUFFER_SIZE_SMALL)
* sambamba (SAMBAMBA_MARKDUP_OPTS)
* samtools (SAMPESORT_MEMSIZE)
* picard (PICARD_MARKDUP_JVM_OPTS)

Other relevant options are

* The resource requirements depend on the workflow variant that is used (e.g. whether biobambam's bamsort or samtools sort is used to sort the BAM file). Use a fast duplication marker. We found sambamba-0.5.9 to be optimal. Use the `markDuplicatesVariant` variable.
* If you have a large and fast local filesystem thes `useRoddyScratchAsBigFileScratch` to true and set `scratchBaseDirectory` in the `applicationProperties.ini` to a path on that filesystem. This will speed up all temporary file IO.
* Mbuffer is used to buffer short timescale throughput fluctuation and for copying data to multiple output (named) pipes.

# BAM and FASTQ Usage

The workflow was originally designed to retrieve the files -- FASTQs and possible BAMs -- from the filesystem by predicting their names using the filename patterns. Usually, the workflow was run on one directory (both input and output directory) and the workflow made sure that all files it expects are present.

It depends on the set of parameters, which BAM file is used as input, e.g. when doing incremental merge.

* If the `bam` parameter is set to some file, an existing target BAM file (derived from the filename pattern in the XML) will be rescued by renaming with a date suffix.
* If `fastq_list` is set (to a semicolon-separated list of files) FASTQ files matching filename patterns are ignored. However, read groups already in the BAM file (no matter whether provided via `bam` or via filename pattern matching) are ignored.
* `bam` and `fastq_list` can be set together. Then all BAM and FASTQ files matching filename patterns are ignored.
* If `useOnlyExistingTargetBam` is set to true, then all existing FASTQ files matching filename patterns are ignored.
* Using `useOnlyExistingTargetBam=true` with 'bam' or 'fastq_list' set is considered a configuration error. The workflow will check this and stop before the job submission.

## Read Group Identifiers

Read group IDs in BAMs are determined (input files) from or stored in (output files) the `ID` attribute in `@RG` header lines. Usually, read group IDs in FASTQ files are determined from filenames using the patterns `${RUN}_${LANE}`. With a metadata input table you can provide FASTQ files with arbitrary file names, because the metadata is taken from the table's columns.

# Running the Workflow

The most important parameters are:

| Parameter  | Example      | Description                            |
|------------|--------------|----------------------------------------|
| INDEX_PREFIX | /path/to/assembly/assembly.fa | The path to the fasta file with the assembled genome(s). Note that the BWA index needs to be this directly and use the string 'assembly.fa' as prefix |
| CHROM_SIZES_FILE | /path/to/assembly/sizes.tsv | A two-column TSV file with chromosome identifiers (1) and number of bases (2). Usually you want the number of bases just be from the set {A, T, C, G}, to ignore in the statictics all lower-quality bases in the genome. | 
| CHR_PREFIX | chrMmu       | With xenograft data used to discern 'matching' and 'nonmatching' identifiers which match the /$CHR_PREFIX$chr$CHR_SUFFIX/ pattern or not, respectively. Also used for WGBS. |
| CHR_SUFFIX | _hsa      | See CHR_PREFIX. |
| CHR_GROUP_NOT_MATCHING | nonmatching | See CHR_PREFIX. |
| CHR_GROUP_MATCHING | matching | See CHR_PREFIX. |
| CHROMOSOME_INDICES | "( 1 2 3 )" | Needed for the WGBS workflow to select chromosomes to be processed. The should be a quoted bash array, i.e. with spaces as element separators and including the parentheses. |

A full description of all options in the different workflows can be found in the XML files in `resources/configurationFiles`. Note that workflow configurations inherit from each other in the order "WGS" <- "WES" <- "WGBS". Thus the WGS configuration (analysisQc.xml) contains variables that are overridden by values in the WES configuration (analysisExome.xml), and so forth.

## Xenograft

The WGS and WES workflows can deal with xenograft data. To process xenograft data you need a combined FASTA file and genome index plus a matching "chromosome sizes file". Thus, `INDEX_PREFIX` and `CHROM_SIZES_FILE` need to be set to the paths of the FASTA file (and BWA index) and the "stats" containing the chromosome sizes. Note that for xenograft the `CHROM_SIZES_FILE` should include also the chromosomes of the host (mouse).

If you want to have species-specific coverage you additionally need to set some variables. The chromosomes of one of the two species needs to be prefixed (`CHR_PREFIX`) and/or suffixed (`CHR_SUFFIX`). For instance you may use a FASTA file with the human chromosomes without (explicitly configured) prefixes, e.g. with chromosomes 1, 2, ..., X, Y, MT and mouse chromosomes prefixed by 'chrMmu'. In this situation use the following configuration values:

* `CHR_PREFIX`=chrMmu
* `CHR_GROUP_MATCHING`=mouse
* `CHR_GROUP_NOT_MATCHING`=human

This will do an alignment of all reads against the combined genomes and additionally collect statistics in the quality-control JSON for the chromosomes in human (without the prefix and suffix) and mouse (matching the "chrMmu" prefix) groups.

*NOTE:* Currently, the WGBS workflow variant uses the `CHR_PREFIX`-variable for another purpose and, therefore, can not collect dedicated statistics for xenograft data.   

# Release Branches

Various versions are or have been in production mode at the DKFZ/ODCF. These often have dedicated release branches named "ReleaseBranch_$major.$minor\[.$patch\]" in which only certain changes have been made:

  * No changes that alter previous output.
  * New important features are sometimes backported -- as long as they do not change the previous results.
  * Bugfixes that allow running the workflow on some data on which it previously crashed, but that do not alter the existing output, are included.
  
Note that ReleaseBranch_1.2.182 is __not the newest branch, but the oldest__! It was derived from a very old version of the workflow (QualityControlWorkflows_1.0.182) at a time where the versioning system was not fixed to [semver 2.0](https://semver.org/).