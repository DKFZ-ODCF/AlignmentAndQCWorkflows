# Alignment and Quality Control Plugin for Roddy

This plugins contains a series of alignment and quality control related workflows:
- PanCancer alignment workflow for whole genome and exome
- Postmerge QC workflow for whole genome and exome
- Bisulfite core workflow

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
conda env create -n AlignmentAndQCWorkflows -f $PATH_TO_PLUGIN_DIRECTORY/resources/qcAnalysis/environments/conda.yml
```

The name of the Conda environment is arbitrary but needs to be consistent with the `condaEnvironmentName` variable. The default for that variable is set in `resources/configurationFiles/qcAnalysis.xml`.

### Disclaimer

The software versions in the Conda environment are not exactly the same as the ones in our local HPC infrastructure at the German Cancer Research Center. Notable differences are:

* biobambam: Instead of 0.0.148 the Conda environment contains 2.0.79. As long as you do not select `markDuplicatesVariant=biobambam` this won't be a problem, as biobambam is only used for sorting BAMs.
* bwa: For the WGBS workflow we currently use a patched version of BWA that does not check for the "/1" and "/2" first and second read marks. This version is not available in BioConda and thus the WGBS workflow won't work with the Conda environment.