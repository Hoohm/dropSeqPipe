This pipeline is dependent on conda.

### Step 1: Download and install miniconda3
First you need to download and install miniconda3:

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


### Step 2: Clone the workflow

Clone the worflow
```
git clone https://github.com/Hoohm/dropSeqPipe.git
```

### Step 3: Install snakemake

```
conda install -c bioconda -c conda-forge snakemake
```
 
Next step is config files completion

[Complete the config.yaml](https://github.com/Hoohm/dropSeqPipe/wiki/Create-config-files) with the missing information

### UPDATES: How to update the pipeline

Go to your experiment folder, then pull.
```
git pull https://github.com/Hoohm/dropSeqPipe.git
```

If you want to update files/plots based on the updates you can use this command:
```
snakemake -R `snakemake --list-codes-changes`
```
This will update all the files that would be modified by the changes in the code (rules or script). Depending on how much and where the changes have been made, this might rerun the whole pipeline.