[![Snakemake](https://img.shields.io/badge/snakemake-≥4.1.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/Hoohm/dropSeqPipe.svg?branch=master)](https://travis-ci.org/Hoohm/dropSeqPipe)

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

## [0.31]
### Changed
- Fixed error for STAR index generation. It crashed saying it couldn't write in folder.
- Fixed a missing plot for plot_knee_plot_whitelist.
- Input files for the STAR_align rule have been changed. Adding samples in an already aligned experiment with a different R2 length, will only align the new data and not realign the old one.
- Split reads and barcodes multiqc reports for qc step.
- Modified a few rules to follow the guidelines for [snakemake workflows](https://github.com/snakemake-workflows/docs)
- Fixed an issue where snakemake would crash on clusters if using `expand()` on fixed variables such as `annotation_prefix`. Now using normal python formatting.
- Changed the config.yaml parameters names to lowercase and hyphens! Software specific variables have their original style making it easier to search in manuals. You will have to either copy the new config.yaml from the templates or modify your own accordingly.
- cell-barcode-edit-distance changed to what it actually is, UMI-edit-distance.
- Updated all the envs to fix bugs.
- Fixed a bug where the mixed species would not run properly.

### Added
- Added ggpubr in environment.yaml file.
- Added a `templates` folder which will hold `config.yaml`, `samples.csv`, `cluster.yaml` as well as adapters files. This will also help cloning the repository without overwritting your own config.yaml file when updating the pipeline.
- Added the possibility of using your own adapters fasta file for trimmomatic. To use it, please refer to the [WIKI](https://github.com/Hoohm/dropSeqPipe/wiki/Create-config-files#filter)
- Added fastqc, multiqc, STAR wrappers. You have now to use the `--use-conda` option to run the pipeline.
- Added cluster recommendations on the wiki.
- Added Localrules for certain rules. This allows to run low ressource rules on the host computer instead of nodes when using clusters.
- genomeChrBinNbits will be calculated automacially for STAR.
- Exposed all variables for trimmomatic in config.yaml under trimming.

### Removed
- png plots have been removed. It was causing some issues on clusters with cairo. Usability is more important than png plots to me.


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
