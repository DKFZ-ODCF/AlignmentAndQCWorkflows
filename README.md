# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/eilslabs/Roddy) workflows:

- PanCancer alignment workflow for whole genome (WGS) and exome (WES)
- Bisulfite core workflow (WGBS) using [methylCtools](https://github.com/hovestadt/methylCtools)

These are basically BWA alignment (bwa mem) workflows with plenty of additional quality control steps. QC-steps acting on the BAM files are mostly done through piping the input data into multiple QC tools simultaneously to achieve a high performance. All workflows can be run for each sample or with combinations of tumor- and control-samples.

> <table><tr><td><a href="https://www.denbi.de/"><img src="docs/images/denbi.png" alt="de.NBI logo" width="300" align="left"></a></td><td><strong>Your opinion matters!</strong> The development of this workflow is supported by the <a href="https://www.denbi.de/">German Network for Bioinformatic Infrastructure (de.NBI)</a>. By completing <a href="https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=AlignmentAndQCWorkflows">this very short (30-60 seconds) survey</a> you support our efforts to improve this tool.</td></tr></table>


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

The current implementation of the WGBS workflow uses [methylCtools](https://github.com/hovestadt/methylCtools) requires a patched BWA version. Note that the methylCtools version in this repository differs from the one in the official Github repository.

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
* `bam` and `fastq_list` can be set together. Then all only the explicitly set files are used but other BAM and FASTQ files matching filename patterns are ignored.
* If `useOnlyExistingTargetBam` is set to true, then all existing FASTQ files matching filename patterns are ignored.
* Using `useOnlyExistingTargetBam=true` with 'bam' or 'fastq_list' set is considered a configuration error. The workflow will check this and stop before the job submission.

## Read Group Identifiers

Read group IDs in BAMs are determined (input files) from or stored in (output files) the `ID` attribute in `@RG` header lines. Usually, read group IDs in FASTQ files are determined from filenames using the patterns `${RUN}_${LANE}`. With a metadata input table you can provide FASTQ files with arbitrary file names, because the metadata is taken from the table's columns.

# Configuration

All configuration variables are documented in the XML files in the `resources/configurationFiles/` directory. There is one XML for each workflow. Note that the workflows depend on each other, i.e. the WES and WGBS workflows are extension of the WGS workflow -- this can be recognized from the `imports` attribute of the top-level configuration tag in the XMLs. This means that most options of the WGS workflow also affect the other two workflows. Conversely, settings in the WGBS and WES workflow may override those in the WGS workflow. Some processing steps, notably those of the ACEseq quality control (QC) are not valid for the WES workflow. Note that the plugin depends on the [COWorkflowBasePlugin](https://github.com/TheRoddyWMS/COWorkflowsBasePlugin), which has its own configurations affecting this alignment plugin.

See the [Roddy](https://github.com/TheRoddyWMS/Roddy) documentation for a description of how to configure and run workflows.

The most important parameters are:

| Parameter  | Example      | Description                            |
|------------|--------------|----------------------------------------|
| INDEX_PREFIX | /path/to/assembly/assembly.fa | The path to the fasta file with the assembled genome(s). Note that the BWA index needs to be this directly and use the string 'assembly.fa' as prefix |
| CHROM_SIZES_FILE | /path/to/assembly/sizes.tsv | A two-column TSV file with chromosome identifiers (1) and number of bases (2). Usually you want the number of bases just be from the set {A, T, C, G}, to ignore all lower-quality bases in the genome in the statistics. | 
| CHR_PREFIX | chrMmu       | With xenograft data used to discern 'matching' and 'nonmatching' identifiers which match the /$CHR_PREFIX$chr$CHR_SUFFIX/ pattern or not, respectively. Also used for WGBS. |
| CHR_SUFFIX | _hsa      | See CHR_PREFIX. |
| CHR_GROUP_NOT_MATCHING | human | See CHR_PREFIX. Default: "nonmatching" |
| CHR_GROUP_MATCHING | mouse | See CHR_PREFIX. Default" "matching" |
| CHROMOSOME_INDICES | "( 1 2 3 )" | Needed for the WGBS workflow to select chromosomes to be processed. This should be a quoted bash array, i.e. with spaces as element separators and including the parentheses. |
| ??? | Trimming is done before alignment. |
| ADAPTOR_TRIMMING_OPTIONS_1 | ? | |
| ADAPTOR_TRIMMING_OPTIONS_1 | ? | |

A full description of all options in the different workflows can be found in the XML files in `resources/configurationFiles`. Note that workflow configurations inherit from each other in the order "WGS" <- "WES" <- "WGBS". Thus the WGS configuration (analysisQc.xml) contains variables that are overridden by values in the WES configuration (analysisExome.xml), and so forth.


## Xenograft

The WGS and WES workflows can deal with xenograft data simply by aligning against the combined genome. Thus to process xenograft data you need a FASTA file, a genome index, and a matching "chromosome sizes file" -- all with both the host's and xenografted species's genomes. Make sure that the chromosomes from both species have different identifiers, e.g. by pre- or suffixing one of sets of the chromosome names, e.g. with chrMmu or whatever is appropriate.

If you want to have species-specific coverage you additionally need to set some variables. The chromosomes of one of the two species need to be prefixed (`CHR_PREFIX`) and/or suffixed (`CHR_SUFFIX`). For instance you may use a FASTA file with the human chromosomes without (explicitly configured) prefixes, e.g. with chromosomes 1, 2, ..., X, Y, MT and mouse chromosomes prefixed by 'chrMmu'. In this situation use the following configuration values:

* `CHR_PREFIX`=chrMmu
* `CHR_GROUP_MATCHING`=mouse
* `CHR_GROUP_NOT_MATCHING`=human

The result will be that in the files `$sample_$pid(_targetExtract)?.rmdup.bam.DepthOfCoverage(_Target)?_Grouped.txt` two lines are created called "matching" (or "mouse" in the example) and "nonmatching" (or "human in the example). Additionally, these values are collected into the `_qualitycontrol.json` files.

*NOTE:* Currently, the WGBS workflow variant uses the `CHR_PREFIX`-variable for another purpose and, therefore, can not collect dedicated statistics for xenograft data.   

## Information for specific protocols

The plugin contains three related workflows for WGS, WES and WGBS data. The way to invoke a specific workflows is to set the `availableAnalyses` section in the project configuration to the configuration name of the desired workflow (i.e. the name of the configuration file in the plugin's `resources/configurationFiles` directory). E.g. the following would define an analysis "WGBS" referring to the bisulfite workflow configuration: 

```xml
<availableAnalyses>
    <analysis id="WGBS" configuration="bisulfiteCoreAnalysis"/>
</availableAnalyses>
```

## Whole Genome Sequencing (WGS) 

The WGS variant does some GC- and replication-timing bias corrections for the coverage estimates, as are described is the documentation of the [ACEseq workflow](https://aceseq.readthedocs.io/en/latest/methods.html#gc-replication-timing-bias-correction).

![WGS job structure](http://www.plantuml.com/plantuml/proxy?cache=no&src=https://raw.github.com/DKFZ-ODCF/AlignmentAndQCWorkflows/master/docs/images/jobs-wgs.puml)

## Whole Exome Sequencing (WES)

For exome sequencing the QC statistics need to account for the target regions only, otherwise the estimates would be widely off any relevant value.

![WES job structure](http://www.plantuml.com/plantuml/proxy?cache=no&src=https://raw.github.com/DKFZ-ODCF/AlignmentAndQCWorkflows/master/docs/images/jobs-wes.puml)

exome-related options



## Whole Genome Bisulfite Sequencing (WGBS)

The WGBS variant does bisulfite calling on the fly with a patched version of [methylCtools](https://github.com/hovestadt/methylCtools) that is included in this repository. The patched version was extended to also be able to handle tagmentation-based data (see next sections).

Here is the overall structure of the WGBS workflow. The additional "library"-merging step is only invoked if there are multiple libraries provided, which is usually only the case for tagmentation data.

![WGBS job structure](http://www.plantuml.com/plantuml/proxy?cache=no&src=https://raw.github.com/DKFZ-ODCF/AlignmentAndQCWorkflows/master/docs/images/jobs-wgbs.puml)

The WGBS workflow is invoked if the "bisulfiteCoreAnalysis" configuration is referenced in the `<availableAnalyses>` section of the config file. Here a complete example for a project configuration for a standard WGBS analysis. 

```xml
<configuration
        configurationType="project"
        name="configurationName"
        description="The description">

    <subconfigurations>
        <configuration name="config" usedresourcessize="xl">
            <availableAnalyses>
                <analysis id="WGBS" configuration="bisulfiteCoreAnalysis"/>
            </availableAnalyses>
        </configuration>
    </subconfigurations>
</configuration>

```

### Swift Biosciences ACCEL-NGS 1S PLUS & METHYL-SEQ

The protocol produces a second read (R2) fragment-end of on average 8 bp containing non-genomic sequences with low complexity. As these sequences are not of genomic origin they should be trimmed off. Because of read through with fragments shorter than the read length, the advise by Swift Biosciences is to trim off this non-genomic 10 bp from *both* fragment ends (compare [here](https://www.westburg.eu/uploads/fckconnector/6f34a0c2-8c61-4ef2-b79d-fe45db09ba51/2888031684)). Like for the WGS workflow, the trimming off is done by trimmomatic before the actual alignment and can be customized as described by adapting the option `ADAPTOR_TRIMMING_OPTIONS_1`. The variable `ADAPTOR_TRIMMING_OPTIONS_1_SwiftAccelNgs` has the correct trimming parameters for the Swift ACCEL-NGS protocol predefined.  

```xml
<cvalue name="IS_TAGMENTATION" value="false"/>
<cvalue name="ADAPTOR_TRIMMING_OPTIONS_1" value="${ADAPTOR_TRIMMING_OPTIONS_1_SwiftAccelNgs}"/>
```

Note that here `IS_TAGMENTATION` is set to false, so no additionally ignore 9 bp during bisulphite calling. 

### WGBS-Tagmentation

WGBS-tagmentation ([Wang _et al._, 2013](https://doi.org/10.1038/nprot.2013.118)) produces about 9 bp of genomic sequences on both fragment ends that show a conversion bias. Because these sequences are genomic, they contain information for the alignment and should not get trimmed off completely. However, because they are biased they need to be ignored during the bisulphite calling. This is the function of the patch of [methylCtools](https://github.com/hovestadt/methylCtools), to ignore the biased 9 bp. Therefore, for tagmentation you need to set 

```xml
<cvalue name="IS_TAGMENTATION" value="true"/>
``` 

in your configuration.

Note that tagmentation data is based on independently amplified libraries, which makes it necessary to do independent duplication marking for each individual library before merging everything into a final merged-BAM.

### PBAT

The Post-Bisulfite Adapter Tagging (PBAT; [Miura _et al._, 2012](https://doi.org/10.1093/nar/gks454)) protocol produces undirectional read pairs. 

By setting `reorderUnidirectionalWGBSReadPairs` the a read-reordering script will be run that decides based on the relative frequencies of TC and AG dinucleotides in both reads, what is the most likely correct orientations of the reads, and may then swap the two reads. Reads that cannot be unambiguously classified are currently dropped. Note that after the swapping, the read-numbers of swapped reads are reversed: What was R1 in the input FASTQ will be R2 in the output BAM, and vice versa. The original script for swapping, including a documentation of the underlying ideas, can be found [here](https://github.com/cimbusch/TWGBS.git).

```xml
<cvalue name="IS_TAGMENTATION" value="false"/>
<cvalue name="reorderUnidirectionalWGBSReadPairs" value="true"/>

```

# Change Logs

See [here](Changelog.md) for the general change logs of the master branch.

## Release Branches

Various versions are or have been in production mode at the DKFZ/ODCF. These often have dedicated release branches named "ReleaseBranch_$major.$minor\[.$patch\]" in which only certain changes have been made:

  * No changes that alter previous output.
  * New important features are sometimes backported -- as long as they do not change the previous results.
  * Bugfixes that allow running the workflow on some data on which it previously crashed, but that do not alter the existing output, are included.
  
> Note that [ReleaseBranch_1.2.182](../../tree/ReleaseBranch_1.0.182) is __not the newest branch, but the oldest__! It was derived from a very old version of the workflow ([QualityControlWorkflows_1.0.182](../../tree/ReleaseBranch_1.0.182)) at a time where the versioning system was not fixed to [semver 2.0](https://semver.org/).

