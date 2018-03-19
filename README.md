[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

Description
------------------
This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It allows to go from raw data of your Single Cell RNA seq experiment until the final count matrix with QC plots along the way.

This is the tool we use in our lab to improve our wetlab protocol as well as provide an easy framework to reproduce and compare different experiments with different parameters.

It uses STAR to map the reads. It is usable for any single cell protocol using two reads where the first one holds the Cell and UMI barcodes and the second read holds the RNA. Here is a non-exhausitve list of compatible protocols:

* Drop-Seq
* SCRB-Seq
* 10x Genomics
* DroNc-seq

This package is trying to be as user friendly as possible. One of the hopes is that non-bioinformatician can make use of it without too much hassle. It will still require some command line execution, this is not going to be fully interactive package.


Latest changes
-----------------

Version 0.31




installation
------------

The installation process can be found in the [Installation](https://github.com/Hoohm/dropSeqPipe/wiki/Installation) section of the wiki.


Future implementations
---------------------------
I'm actively seeking help to implement the points listed bellow. Don't hesitate to contact me if you wish to contribute.

* Create a sharing platform where quality plots/logs can be discussed and troubleshooted.
* Create a full html report for the whole pipeline
* Multiqc module for drop-seq-tools
* Conda package for drop-seq-tools
* RData object of all the summary data and plots so that you can create your own report.
* Implement an elegant "preview" mode where the pipeline would only run on a couple of millions of reads and allow you to have an approximated view before running all of the data. This would dramatically reduce the time needed to get an idea of what filters whould be used.

I hope it can help you out in your single cell experiments!

Feel free to comment and point out potential improvements.