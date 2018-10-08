# Config file and sample file
---------------------------

In order to run the pipeline you will need to complete the config.yaml file and the samples.csv file. Both are located in the `templates` folder and should be moved to the root folder of the experiment.

## 1. config.yaml - Executables, system and experiment parameters
The config.yaml contains the paths to the the Drop-Seq tools wrapper, paths to references and parameters for each step of the pipeline.
```
    :*emp-directory: /path/to/temp/or/scratch/folder
    r*q-wrapper: /path/to/drop-seq-tools-wrapper.sh
memory: 4g
META:
    *cies:
    *   - SPECIES_ONE
    *   - SPECIES_TWO
    ratio: 0.2
    reference-file: reference.fasta
    annotation-file: annotation.gtf
    reference-directory: /path/references/files/
FILTER:
    5-prime-smart-adapter: DEPENDS ON THE PROTOCOL
    cell-barcode:
        start:
        end:
        min-quality:
        num-below-quality:
    UMI-barcode:
        start:
        end:
        min-quality:
        num-below-quality:
    trimmomatic:
        adapters-file: /path/to/adapters.fa
        LEADING: 3
        TRAILING: 3
        SLIDINGWINDOW:
            windowSize: 4
            requiredQuality: 20
        MINLEN: 20
        ILLUMINACLIP:
            seedMismatches: 2
            palindromeClipThreshold: 30
            simpleClipThreshold: 10
MAPPING:
    STAR:
        outFilterMismatchNmax: 10
        outFilterMismatchNoverLmax: 0.3
        outFilterMismatchNoverReadLmax: 1
        outFilterMatchNmin: 0
        outFilterMatchNminOverLread: 0.66
        outFilterScoreMinOverLread: 0.66
EXTRACTION:
    UMI-edit-distance:
    minimum-counts-per-UMI:


```
Please note the "space" after the colon, is needed for the yaml to work.

## Subsections
### [LOCAL]
* `temp-directory` is the temp or scratch folder with enough space to keep temporary files.
* `dropseq-wrapper` is the wrapper that will call all the drop-seq-tools executables.
* `memory` is the maximum memory allocation pool for a Java Virtual Machine.

### [META]
* `species` is the species used in the experiment. [Required for mixed species experiment]
Ex:
```
 species:
    - MOUSE
    - HUMAN
```

* `ratio` is how much "contamination" from another species you allow to validate them as a species or mixed. 0.2 means you allow a maximum of 20% mixing.
Note: Those species name must reflect names used in the genome.fasta

* `reference-file` is the reference file of your genome (fasta).
* `annotation-file` is the gene annotation file (GTF).
* `reference-folder` is the folder where both of these files are gonna be stored.

### [FILTER]

* `5-prime-smart-adapter` is the 5" smart adapter used in your protocol.
* `cell-barcode and UMI-barcode`: Is the section for cell/umi barcode filtering.
    * `start` is the first base position of your cell/umi barcode.
    * `end` is the last base position of your cell/umi barcode.
    * `min-quality` is the minimum quality filtering for bases in your cell barcodes or UMI to discard reads.
    * `num-below-quality` is the maximum bases under `min_quality`. A value of `0` means all reads that have more than 0 (1 or more) bases under `min-quality` in the cell/umi barcode are discarded.
    
* `trimmomatic`: Is the section for trimming.
* `adapters-file` is the file containing your list of adapters as fasta. you can choose between 6 files in the `templates` folder, add any sequence to existing files or provide your own custom one.
    * NexteraPE-PE.fa 
    * TruSeq2-PE.fa
    * TruSeq2-SE.fa
    * TruSeq3-PE-2.fa
    * TruSeq3-PE.fa
    * TruSeq3-SE.fa
Provide the path to the file you want to use for trimming. If you want to add custom sequences or create a complete new one, I would advise to store it in the ROOT folder of the experiment. This will ensure that your custom file will not be overwritten if you update the pipeline.

Example: `NexteraPE-PE.fa`

* `LEADING` (default:3) quality: Specifies the minimum quality required to keep a base.
* `TRAILING` (default:3) quality: Specifies the minimum quality required to keep a base.
* `SLIDINGWINDOW`: 
    * `windowSize` (default:4) specifies the number of bases to average across.
    * `requiredQuality` (default:20) specifies the average quality required.
* `MINLEN` (default:20) Specifies the minimum length of reads to be kept.
* `ILLUMINACLIP`:
    * `seedMismatches` (default:2) specifies the maximum mismatch count which will still allow a full match to be performed.
    * `palindromeClipThreshold` (default:30) specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
    * `simpleClipThreshold` (default:10) specifies how accurate the match between any adapter etc. sequence must be against a read.
All of the values for trimmomatic are the default ones. For details about trimmomatic parameters and what they do, please refer to the [trimmomatic web page](http://www.usadellab.org/cms/?page=trimmomatic).


### [MAPPING]
* `STAR`
    * `outFilterMismatchNmax` (default:10) is the maximum number of mismatches allowed.
    * `outFilterMismatchNoverLmax` (default:0.3) is the maximum ratio of mismatched bases that mapped.
    * `outFilterMismatchNoverReadLmax` (default:1.0) is the maximum ratio of mismatched bases of the whole read.
    * `outFilterMatchNmin` (default:0) is the minimum number of matched bases.
    * `outFilterMatchNminOverLread` (default:0.66) alignment will be output only if the ratio of matched bases is higher than or equal to this value.
    * `outFilterScoreMinOverLread` (default:0.66) alignment will be output only if its ratio score is higher than or equal to this value.

All of the values for STAR are the default ones. For details about STAR parameters and what they do, please refer to the [STAR manual on git](https://github.com/alexdobin/STAR/tree/master/doc).


### [EXTRACTION]
* `UMI-edit-distance` This is the maximum manhattan distance between two UMI barcode when extracting count matrices.
* `min-count-per-umi` is the minimum UMI/Gene pair needed to be counted as one.


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

Next step [downloading the reference files](https://github.com/Hoohm/dropSeqPipe/wiki/Reference-Files)