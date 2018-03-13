[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

Description
------------------
This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It allows to go from raw data of your dropSeq/scrbSeq experiment until the final count matrix with QC plots along the way.
This is the tool we use in our lab to improve our wetlab protocol as well as provide an easy framework to reproduce and compare different experiments with different parameters.

It uses STAR to map the reads. Is is working for single cell UMI RNAseq data such as:

* Drop-Seq
* SCRB-Seq
* 10x Genomics
* DroNc-seq

This package is trying to be as user friendly as possible. One of the hopes is that non-bioinformatician can make use of it without too much hassle. It will still require some command line execution, this is not going to be an interactive package.


Latest changes
-----------------

Version 0.3 is out. The python package got changed to a snakemake-workflow. Main differences listed bellow

* It now uses [conda](https://conda.io/docs/) for packages/installation. This makes the installation process easier and will allow to use the pipeline in any cluster even if you don't have admin rights.
* All of the [snakemake](http://snakemake.readthedocs.io/en/latest/) features are now fully exposed whereas before they had to be exposed manually before.
* [Multiqc](http://multiqc.info/) is now properly integrated in the pipeline.
* Local, software configurations and sample specific informations are now separated in two different files, `config.yaml` and `samples.csv`



installation
------------

The installation process can be found in the [Installation](https://github.com/Hoohm/dropSeqPipe/wiki/Installation) section of the wiki.


Future implementations
---------------------------
I'm actively seeking help to implement the points listed bellow. Don't hesitate to contact me if you wish to contribute.

* Multiqc module for drop-seq-tools
* Conda package for drop-seq-tools
* Multithreading for java tasks
* RData object of all the summary data and plots so that you can create your own report.
* Implement an elegant "preview" mode where the pipeline would only run on a couple of millions of reads and allow you to have an approximated view before running all of the data.

I hope it can help you out in your drop-seq experiment!

Feel free to comment and point out potential improvements.