# Config file and sample file
---------------------------

In order to run the pipeline you will need to complete the config.yaml file and the samples.csv file. Both are located in the `templates` folder , should be moved to the root folder of the experiment and filled in for missing entries before running the pipeline.

The goal for this is to provide the config.yaml when you finally upload the data to a repository for a publication as well as the pipeline version. This provides other users to ability to rerun the processing from scratch exactly as you did. This is possible because snakemake will download and create the exact same environnment for each rule using the envs files provided with the pipeline.

## 1. config.yaml - Executables, system and experiment parameters
The config.yaml contains all the necessary parameters and paths for the pipeline.
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
Please note the "space" after the colon, is needed for the yaml to work.

## Subsections

### [CONTACT]
* `email` and `person` This is not requested. You can provide the e-mail and name address of the person who processed the data using this configuration. Ideally you should provide the config.yaml with the data repository to allow people to rerun the data using dropSeqPipe.

### [LOCAL]
* `temp-directory` is the temp or scratch folder with enough space to keep temporary files.
* `memory` is the maximum memory allocation pool for a Java Virtual Machine.
* `raw_data` is the folder containing all your raw fastq.gz files.
* `results` is the folder that will contain all the results of the pipeline.

### [META]
* `species` is where you list the species of your samples. It can be a mixed experiment with two entries. 
* `SPECIES_ONE` can be for example: mus_musculus, homo_sapiens, etc... It has to be the name used on ensembl for automatic download to work.
* `build` is the genome build number.
* `release` is the annotation release number.
* `SPECIES_TWO` can be your second species.
* `ratio` is how much "contamination" from another species you allow to validate them as a species or mixed. 0.2 means you allow a maximum of 20% mixing.
* `reference-directory` is where you want to store your references files.
* `gtf_biotypes` is the gtf_biotypes.yaml file containing the selection of biotypes you want to keep for your gene to read attribution. Using less biotypes may decrease your multimapping counts.

### [FILTER]
* `barcode_whitelist` is the filename of your whitelist fi you have one. Well plate base protocols often have one.
* `5-prime-smart-adapter` is the 5" smart adapter used in your protocol.
* `cell-barcode and UMI-barcode`: Is the section for cell/umi barcode filtering.
    * `start` is the first base position of your cell/umi barcode.
    * `end` is the last base position of your cell/umi barcode.    
* `cutadapt`: Is the section for trimming.
* `adapters-file` is the file containing your list of adapters as fasta. you can choose between 6 files in the `templates` folder, add any sequence to existing files or provide your own custom one.
    * NexteraPE-PE.fa 
    * TruSeq2-PE.fa
    * TruSeq2-SE.fa
    * TruSeq3-PE-2.fa
    * TruSeq3-PE.fa
    * TruSeq3-SE.fa
Provide the path to the file you want to use for trimming. If you want to add custom sequences or create a complete new one, I would advise to store it in the ROOT folder of the experiment. This will ensure that your custom file will not be overwritten if you update the pipeline.

Example: `NexteraPE-PE.fa`
* `R1` lists the options for read1 (cell barcode and umi) filtering/trimming
    * `quality-filter` is the minimum mean score of the sliding window for quality filtering.
    * `maximum-Ns` how many Ns you allow in the cell barcode and umi barcode. By default it is one because we want to be able to collapse barcodes that have one mismatch.
    * `extra-params` if you usually add extra paramters to cutadapt, you can do it here. *Only for experienced cutadapt users*.
* `R2` lists the options for read2 (mRNA) filtering/trimming
that have one mismatch.
    * `maximum-length` is the maximum length of your mRNA read before alignement.
    * `extra-params` if you usually add extra paramters to cutadapt, you can do it here. *Only for experienced cutadapt users*.
For more information about trimming and filtering please visit the [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) website.

### [MAPPING]
* `STAR`
    * `genomeChrBinNbits` is a value used for index generation in STAR. The formula is min(18,int(log2(genomeLength/referenceNumber)))
    * `outFilterMismatchNmax` (default:10) is the maximum number of mismatches allowed.
    * `outFilterMismatchNoverLmax` (default:0.3) is the maximum ratio of mismatched bases that mapped.
    * `outFilterMismatchNoverReadLmax` (default:1.0) is the maximum ratio of mismatched bases of the whole read.
    * `outFilterMatchNmin` (default:0) is the minimum number of matched bases.
    * `outFilterMatchNminOverLread` (default:0.66) alignment will be output only if the ratio of matched bases is higher than or equal to this value.
    * `outFilterScoreMinOverLread` (default:0.66) alignment will be output only if its ratio score is higher than or equal to this value.

All of the values for STAR are the default ones. For details about STAR parameters and what they do, please refer to the [STAR manual on git](https://github.com/alexdobin/STAR/tree/master/doc).

### [EXTRACTION]
* `LOCUS` are the overlapping regions that reads overlap and are counted in the final expression matrix. Possible values are `CODING`, `UTR`, `INTRON`
* `UMI-edit-distance` This is the maximum manhattan distance between two UMI barcode when extracting count matrices.
* `min-count-per-umi` is the minimum UMI/Gene pair needed to be counted as one.
* `strand-strategy` `SENSE` defines that you only count genes where the forward strand mapped to the forward region on the DNA. Other possibilities are `ANTISENSE` (only count reads that mapped on the opposite strand) or `BOTH` (count all).

# 2. samples.csv - Samples parameters
This file holds the sample names, expected cell numbers and read length for each sample.
The file has to have this format:

```
samples,expected_cells,read_lengths,batch
sample_name1,500,100,Batch1
sample_name2,500,100,Batch2
```

* `expected_cells` is the amount of cells you expect from your sample.
* `read_length` is the read length of the mRNA (Read2). This is necessary for STAR index generation
* `batch` is the batch of your sample. If you are added new samples to the same experiment, this is typically a good place to add the main batch.

`Note:` You can add any other column you wish here, it won't affect the pipeline and you can use it later on in your analysis.

Finally, you can now [run the pipeline](https://github.com/Hoohm/dropSeqPipe/wiki/Running-dropSeqPipe)

or

Create a [custom reference](https://github.com/Hoohm/dropSeqPipe/wiki/Reference-Files)

