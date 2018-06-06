# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/eilslabs/Roddy) workflows:

- PanCancer alignment workflow for whole genome (WGS) and exome (WES)
- Bisulfite core workflow (WGBS)

These are basically BWA alignment (bwa mem) workflows with plenty of additional quality control steps. QC-steps acting on the BAM files are mostly done through piping the input data into multiple QC tools simultaneously to achieve a high performance. All workflows can be run for each sample or with combinations of tumor- and control-samples.

# Configuration

All configuration variables are documented in the XML files in the `resources/configurationFiles/` directory. There is one XML for each workflow. Note that the workflows depend on each other, i.e. the WES and WGBS workflows are extension of the WGS workflow -- this can be recognized from the `imports` attribute of the top-level configuration tag in the XMLs. This means that most options of the WGS workflow also affect the other two workflows. Furthermore, settings in the WGBS and WES workflow may override those in the WGS workflow. Some processing steps, notably those of the ACEseq quality control (QC) are not valid for the WES workflow. 

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

