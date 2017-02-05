This small pipeline allows you to run the basic steps to align and extract expression from a drop-seq experiment.

Before using it you will need some programs:

1. [Snakemake](https://snakemake.readthedocs.io/en/latest/) (based on python)
2. [R](https://cran.r-project.org/)
3. [STAR aligner](https://github.com/alexdobin/STAR)
4. [Drop-seq tools (1.12)](http://mccarrolllab.com/dropseq/)
5. [Picard tools](https://broadinstitute.github.io/picard/)
6. [jsonlite R package](https://cran.r-project.org/web/packages/jsonlite/index.html)

Before running the pipeline you will also need a reference genome as well as the GTF (needed for the ReFlat) and the ReFlat. This is not explained here but you can get the info [here](http://mccarrolllab.com/dropseq/).

In order to run the pipeline you will need to write two config json files.
One in the root folder of the pipeline which will contain the paths to your executables as well as the Drop-Seq tools.
```
{
    "TMPDIR":"/path/to/temp",
    "PICARD":"/path/to/picard/dist/picard.jar",
    "DROPSEQ":"/path/to/Drop-seq_tools-1.12",
    "STAREXEC":"/path/to/STAR/bin/Linux_x86_64/STAR",
    "CORES": X
}
```

I had some issues because I had not enough space on / so I added a temp folder to fix that. Note: If you have the same problem, you have to manually edit all the *.sh files in the drop-seq tools to use this  TMPDIR variable.
* TMPDIR is the temp folder on the disk with enough space
* PICARD is the path to the picard.jar
* DROPSEQ is the path to the folder of Drop-Seq tools
* STAREXEC is the path to the STAR executable
* CORES is the number of cores you want to use in the pipeline (snakemake is great at balancing tasks!)

The other json file should be in the folder containing all your fastq files and should look like that.
```
{
    "Samples": {
        "Sample1":"N701",
        "Sample2":"N702",
        "Sample3":"N703"
        },
    "Primers":{
        "N701":"TCGCCTTA",
        "N702":"CTAGTACG",
        "N703":"TTCTGCCT",
        "N704":"GCTCAGGA"
    },
    "Barcodes":200,
    "GENOMEREF": "/path/to/reference.fa",
    "REFFLAT": "/path/to/reference.refFlat",
    "METAREF": "/path/to/STAR_REF"
}

```

Samples contains a list of the names of your samples. In this example the samples in the folder should look like:
* Sample1_R1.fastq.gz
* Sample1_R2.fastq.gz
* Sample2_R1.fastq.gz
* Sample2_R2.fastq.gz
* Sample3_R1.fastq.gz
* Sample3_R2.fastq.gz

Primers are the common primer used in [the nextera kit](http://seq.liai.org/204-2/).
* Barcodes should be double the amount of cells you expect from your experiment.
* GENOMEREF is the reference fasta of your genome.
* REFLAT is the reference refFlat file needed the pipeline. You can check how to create it in the [Drop-Seq alignement cookbook](http://mccarrolllab.com/dropseq/).
* METAREF is the folder of the STAR index

Once everything is in place, you can run the pipeline using the following command:

`python3 dropseq.py /path/to/your/samples/`


This will create necessary folders in the sample folder and run the pipeline until the knee plot.
Note: The reason why I run the script in three parts is because of the way STAR handles the loading of the reference genome.
The main idea is that I want to load the reference once, process all the samples and then unload the reference.

Future implementations:
* Automated Species plot or Barnyard plot
* Automated extraction of the raw count matrix
* Cluster version (One of the reasons it's based on snakemake)


I hope it can help you out in your drop-seq experiment!

Feel free to comment and point out potential improvements.