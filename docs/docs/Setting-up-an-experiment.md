# Setting up an experiment

The pipeline needs to be installed as shown previously.
In order to run the pipeline, a separate folder (working directory) needs to be created containing raw sequencing data as well as all required configuration files.
dSP will then be run individually on each of the `workingdir` which contain the sample files.
The `workingdir` folder can contain multiple runs (batches) allowing new samples to be added and analysed easily.
Adding new samples, will only run the pipeline on the newly added data and recreate reports as well as plots containing all the samples.
This automatically done by snakemake and will save unecesarry computation performed already.

Example: If running 2 biological conditions with 2 replicates. This makes up for 4 samples. Assuming a simple dropseq protocol with only human cells.
1. Sequence the data, download the respective 8 data files (2 files per sample) and install the pipeline.
2. Run the pipeline with the command: `snakemake --use-conda --cores N --directory WORKING_DIR`. `N` being the number of cores available and `WORKING_DIR` being the folder containing the `config.yaml`, `samples.csv`, adapter file and `gtf_biotypes.yaml`.
3. If new samples have to be added, create a new set of libraries and sequence them (same 2x2 design)
5. Add the new files in the data folder of `WORKING_DIR` and edit the samples.csv to add missing samples.
6. Run the pipeline like the first time `snakemake --use-conda --cores N --directory WORKING_DIR`
7. This will run the new samples only and recreate the reports as well as the yield plots.
8. It is now easy to compare the impact of the changes in the protocol

## Working dir folder preparation

The raw data from the sequencer are stored in the `RAW_DATA` folder of the `WORKING_DIR`, which should be structured as shown below:

```
/path/to/your/WORKING_DIR/
| -- RAW_DATA/
| -- -- sample1_R1.fastq.gz
| -- -- sample1_R2.fastq.gz
| -- -- sample2_R1.fastq.gz
| -- -- sample2_R2.fastq.gz
| samples.csv
| config.yaml
| gtf_biotypes.yaml
|
| custom_adapters.fa
```
Further to the `RAW_DATA` folder the `WORKING_DIR` should contain the `samples.csv`, `config.yaml` , `gtf_biotypes.yaml` and `custom_adapters.fa` which are explained in more detail in the following chapters.


*Note: DropSeq or ScrbSeq contain two files per sample as they are methods using paired end sequencing.
R1 will hold the information of the barcode and UMI, R2 will hold the 3' end of the captured mRNA.*

Optional: If a barcode whitelist is available like in methods such as ScrbSeq add ` barcodes.csv` into the working directory.



## Config file and sample file

In order to run the pipeline complete the `config.yaml` file and the `samples.csv` file.
Templates for both are located in the `templates` folder, which can be copied to the `WORKING_DIR` of the experiment and filled in for missing entries before running the pipeline.

* Note: Ideally the pipeline version and the config files, specifically `config.yaml` and `samples.csv` should be provided with the uploaded data to a repository for a publication.
This provides other users to ability to rerun the processing from scratch with exactly the same parameters.
This is possible because `snakemake` will download and create the exact same conda software environment (if `--use-conda` parameter is used) for each rule using the `envs` files (`envs` files are conda software environment files provided with the pipeline i.e. yaml files contained in the `envs` directory).

### 1. config.yaml
The config.yaml contains all the necessary parameters and paths for the pipeline.

Example

```
CONTACT:
  email: user.name@provider.com
  person: John Doe
LOCAL:
    temp-directory: /tmp
    memory: 4g
    raw_data:
    results:
META:
    species:
        mus_musculus:
            build: 38
            release: 94
        homo_sapiens:
            build: 38
            release: 91
    ratio: 0.2
    reference-directory: /path/to/references/
    gtf_biotypes: gtf_biotypes.yaml
FILTER:
    barcode_whitelist: ''
    5-prime-smart-adapter: AAAAAAAAAAA
    cell-barcode:
        start: 1
        end: 12
    UMI-barcode:
        start: 13
        end: 20
    cutadapt:
        adapters-file: 'adapters.fa'
        R1:
            quality-filter: 20
            maximum-Ns: 0
            extra-params: ''
        R2:
            quality-filter: 20
            minimum-adapters-overlap: 6
            minimum-length: 15
            extra-params: ''
MAPPING:
    STAR:
        genomeChrBinNbits: 18
        outFilterMismatchNmax: 10
        outFilterMismatchNoverLmax: 0.3
        outFilterMismatchNoverReadLmax: 1
        outFilterMatchNmin: 0
        outFilterMatchNminOverLread: 0.66
        outFilterScoreMinOverLread: 0.66
EXTRACTION:
    LOCUS:
        - CODING
        - UTR
    strand-strategy: SENSE
    UMI-edit-distance: 1
    minimum-counts-per-UMI: 0
```
 * Note: The "space" after the colon, is needed for `config.yaml` to work.

#### Subsections

[CONTACT]
* `email` and `person`: Provide the name and e-mail address of the person who processed the data using this configuration, this however is not mandatory.

[LOCAL]
* `temp-directory` is the temp or scratch folder with enough space to keep temporary files.
* `memory` is the maximum memory allocation pool for a Java Virtual Machine.
* `raw_data` is the folder containing all the raw fastq.gz files.
* `results` is the folder that will contain all the results of the pipeline.

[META]
* `species` is where the list the species reference genomes of the samples. It can be one species or two for a mixed species experiment.
* `SPECIES_ONE` can be for example: mus_musculus, homo_sapiens, etc... It has to be the name used on ensembl for automatic download to work.
* `build` is the genome build number.
* `release` is the annotation release number.
* `SPECIES_TWO` can be the second species.
* `ratio` is how much "contamination" from another species is allowed to validate them as a species or mixed. 0.2 means a maximum of 20% mixing is allowed.
* `reference-directory` is where the references files are stored.
* `gtf_biotypes` is the `gtf_biotypes.yaml` file containing the selection of biotypes for the gene to read attribution. Using less biotypes may decrease the multimapping counts.

[FILTER]
* `barcode_whitelist` is the filename of the barcode whitelist.
* `5-prime-smart-adapter` is the 5-prime smart adapter used in the protocol.
* `cell-barcode and UMI-barcode`: Is the section for cell/UMI barcode filtering.
    * `start` is the first base position of the cell/UMI barcode.
    * `end` is the last base position of the cell/UMI barcode.
* `cutadapt`: Is the section for trimming.
* `adapters-file` path to the file containing a list of adapters as fasta  to use for trimming. See section below for further information.

Example: `NexteraPE-PE.fa`
* `R1` lists the options for read1 (cell barcode and UMI) filtering/trimming
    * `quality-filter` is the minimum mean score of the sliding window for quality filtering.
    * `maximum-Ns` how many Ns (N stands for unknown base - instead of A,T,C or G) are allowed in the cell barcode and UMI barcode. By default it is one because to be able to collapse barcodes that have one mismatch.
    * `extra-params` extra parameters to cutadapt can be added here *Only for experienced cutadapt users*.
* `R2` lists the options for read2 (mRNA) filtering/trimming
that have one mismatch.
    * `maximum-length` is the maximum length of the mRNA read before alignment.
    * `extra-params` extra parameters to cutadapt can be added here *Only for experienced cutadapt users*.
For more information about trimming and filtering please visit the [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) website.

[MAPPING]
* `STAR`
    * `genomeChrBinNbits` is a value used for index generation in STAR. The formula is min(18,int(log2(genomeLength/referenceNumber)))
    * `outFilterMismatchNmax` (default:10) is the maximum number of mismatches allowed.
    * `outFilterMismatchNoverLmax` (default:0.3) is the maximum ratio of mismatched bases that mapped.
    * `outFilterMismatchNoverReadLmax` (default:1.0) is the maximum ratio of mismatched bases of the whole read.
    * `outFilterMatchNmin` (default:0) is the minimum number of matched bases.
    * `outFilterMatchNminOverLread` (default:0.66) alignment will be output only if the ratio of matched bases is higher than or equal to this value.
    * `outFilterScoreMinOverLread` (default:0.66) alignment will be output only if its ratio score is higher than or equal to this value.

All of the values for STAR are the default ones. For details about STAR parameters and what they do, please refer to the [STAR manual on git](https://github.com/alexdobin/STAR/tree/master/doc).

[EXTRACTION]
* `LOCUS` are the regions where reads overlap and are counted in the final expression matrix. Possible values are `CODING`, `UTR`, `INTRON`
* `UMI-edit-distance` This is the maximum manhattan distance between two UMI barcode when extracting count matrices.
* `min-count-per-umi` is the minimum UMI/Gene pair needed to be counted as one.
* `strand-strategy` Given the factor `SENSE` this parameter defines that only genes are counted where the forward strand mapped to the forward region on the DNA. Other possibilities are `ANTISENSE` (only counts reads that mapped on the opposite strand) or `BOTH` (count all).

### 2. samples.csv - Samples parameters
This file holds the sample names, expected cell numbers and read length for each sample.
The file has to have this format:

```
samples,expected_cells,read_lengths,batch
sample_name1,500,100,Batch1
sample_name2,500,100,Batch2
```

* `expected_cells` is the number of cells expected from a sample.
* `read_length` is the read length of the mRNA (Read2). This is necessary for STAR index generation
* `batch` Can be used to group samples into batches. If new samples are added to the same experiment, a new batch should be assigned to those samples.

`Note:` Any other column can be added here, it won't affect the pipeline and it can be used later on in the analysis.

### 3. gtf_biotypes.yaml

Here you can select which biotypes you would like to keep in your annotation file. Simply delete the ones you don't want from the file for your experiment and move the folder to the root dir of your experiment.

### 4. custom_adapters.fa

File containing the list of adapters as fasta. The files can be chosen in the `templates` folder and any sequence can be added to existing files or custom files can be provided.
    * NexteraPE-PE.fa
    * TruSeq2-PE.fa
    * TruSeq2-SE.fa
    * TruSeq3-PE-2.fa
    * TruSeq3-PE.fa
    * TruSeq3-SE.fa

If custom sequences are added or new files are created, ist recommended to to store them in the `WORKING_DIR` folder of the experiment. This will ensure that the custom file will not be overwritten if the pipeline is updated.

### 5. Barcode whitelist [optional]

In protocols such as SCRBseq, the expected barcodes sequences are known beforhand.
This pipeline does allow the use of known barcodes instead of a number of expected cells.
In order to use this functionality add a whitelist barcode file and provide the name of the file in the configuration in the section:

```
FILTER:
	barcode_whitelist: name_of_your_whitelist_file
```
The file should be in the WORKING_DIR. Run the pipeline as usual.


### 6. cluster.yaml - Running on clusters [optional]

There is a file in the `templates` called `cluster.yaml`. 
This can be used to modify resources needed for the data. 
Its recommended to move the file to the root of the folder so that it doesn't get replaced by updates.

Below is an example of running on a cluster using the template file `cluster.yaml` for the SLURM Workload Manager used by many of the world's supercomputers and computer clusters (https://slurm.schedmd.com).

```
snakemake --cluster 'sbatch -n {cluster.n}  -t {cluster.time} --clusters=CLUSTERNAME --output={cluster.output}' --jobs N --cluster-config cluster.yaml --use-conda --local-cores C
```

* N: is the number of jobs that are allowed to run at the same time
* C: is the local core of the host machine. A few simple rules are going to be run locally (not sent to nodes) because they are not that heavy (mostly plotting)
* CLUSTERNAME: the name of the cluster that is used.

Note: The default path for cluster logs in the `cluster.yaml` is `logs/cluster/`. If that folder doesn't exist, the cluster can't write and will crash without an error message.

## Choosing the mapping reference

The mapping reference is usually a genome where the reads are aligned to same as for bulk RNA-Seq.
A list of available genomes can be found on ENSEMBL (e.g. release 98):
http://ftp.ensembl.org/pub/release-98/gtf/

From version 0.4 on, reference files are automatically downloaded by the pipeline into a directory specified in the `config.yaml` (`reference-directory: path`).
Mixed references are also downloaded and merged automatically based on the information in the `META` section of `config.yaml` (i.e. if two organism rather than one are listed in the `species: ` section, it will assume it is mixed species).

To use a custom reference and bypass the reference download, create a custom `genome.fa` and `annotation.gtf` file.
Snakemake generates file based on paths.
If a custom reference is used, name it as shown below to enable snakemake to find it.

Here is an example:

Let's assume this is the configuration for the `META` section of `config.yaml`:
```
META:
    species:
        funky_species_name:
            build: A
            release: 1
    ratio: 0.2
    reference-directory: /absolute/path/to/references
    gtf_biotypes: gtf_biotypes.yaml
```

Provide the following files

```
/absolute/path/to/references/funky_species_name_A_1/genome.fa
/absolute/path/to/references/funky_species_name_A_1/annotation.gtf
```

This will stop dropSeqPipe from downloading a new reference.


Once the pipeline has run completely, the `reference-directory` will look like this:

```
genome.fa
annotation.gtf
annotation.refFlat
annotation_reduced.gtf
genome.consensus_introns.intervals
genome.dict
genome.exons.intervals
genome.genes.intervals
genome.intergenic.intervals
genome.rRNA.intervals
STAR_INDEX/SA_read_length/
```

Note: The STAR index will be built based on the read length of the mRNA read (Read2).
If there are different read lengths, it will produce multiple indexes.
