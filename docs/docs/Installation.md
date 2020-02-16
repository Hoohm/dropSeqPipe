# Installation

dropSeqPipe is dependent upon the software management system `conda`.
First step therefore is to install `miniconda3` a small version of `Anaconda` that includes conda, Python, the packages they depend on and a small number of other useful packages such as pip and zlib.


## Step 1: Download and install miniconda3
First download and install miniconda3 by executing the commands below on the command line.

for linux
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

for mac os
```
curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

The command `conda` should now be available and can be tested by running `conda` on the command line.
If the command doesn’t work, more information can be found [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Step 2: Install snakemake

To run the pipeline [snakemake](https://snakemake.readthedocs.io/en/stable/) is needed which can be installed using conda

```
conda install -c bioconda -c conda-forge snakemake
```
The command `snakemake` should now be available and can be tested by running `snakemake` on the command line. If the command doesn’t work, more detailed instruction can be found [here](https://snakemake.readthedocs.io/en/stable/)

## Step 3: Obtain pipeline

In order to download the actual pipeline either download it directly or use [git]( https://www.git-scm.com/) to clone the workflow.

The pipeline can be directly downloaded [here](https://github.com/Hoohm/dropSeqPipe/archive/master.zip)

It is however recommended to download the pipeline by cloning using git rather than downloading the zip archive.
If `git` hasn’t been installed yet `conda` can be used to do so.

```
conda install -c bioconda -c conda-forge git
```

Next the workflow can be obtained using the following command.

```
git clone https://github.com/Hoohm/dropSeqPipe.git
```

This should result in the set-up of a directory called `dropSeqPipe` which contains the pipeline.



## How to update the pipeline

Go to the experiment folder, then pull.

```
git pull
```

To update files/plots, use this command:

```
snakemake -R `snakemake --list-codes-changes`
```
This will update all the files that would be modified by the changes in the code (rules or script). Depending on how much and where the changes have been made, this might rerun the whole pipeline.
