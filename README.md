[![Snakemake](https://img.shields.io/badge/snakemake-≥4.1.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/Hoohm/dropSeqPipe.svg?branch=master)](https://travis-ci.org/Hoohm/dropSeqPipe)

Description
------------------
This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It allows to go from raw data of your Single Cell RNA seq experiment until the final count matrix with QC plots along the way.

This is the tool we use in our lab to improve our wetlab protocol as well as provide an easy framework to reproduce and compare different experiments with different parameters.

It uses STAR to map the reads. It is usable for any single cell protocol using two reads where the first one holds the Cell and UMI barcodes and the second read holds the RNA. Here is a non-exhausitve list of compatible protocols/brands:

* Drop-Seq
* SCRB-Seq
* 10x Genomics
* DroNc-seq
* Dolomite Bio ([Nadia Instrument](https://www.dolomite-bio.com/product/nadia-instrument/))

This package is trying to be as user friendly as possible. One of the hopes is that non-bioinformatician can make use of it without too much hassle. It will still require some command line execution, this is not going to be fully interactive package.


## Authors

* Patrick Roelli ([@Hoohm)](https://github.com/Hoohm))
* Sebastian Mueller ([@seb-mueller)](https://github.com/seb-mueller))
* Charles Girardot ([@cgirardot)](https://github.com/cgirardot))

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/Hoohm/dropSeqPipe/releases).
If you intend to modify and further develop this workflow, fork this reposity. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, once available, its DOI.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml` and the  `samples.tsv` following those [instructions](https://github.com/Hoohm/dropSeqPipe/wiki/Create-config-files)

### Step 3: Execute workflow

All you need to execute this workflow is to install Snakemake via the [Conda package manager](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda). Software needed by this workflow is automatically deployed into isolated environments by Snakemake.

Test your configuration by performing a dry-run via

    snakemake --use-conda -n --directory $WORKING_DIR

Execute the workflow locally via

    snakemake --use-conda --cores $N --directory $WORKING_DIR

using `$N` cores on the `$WORKING_DIR`. Alternatively, it can be run in cluster or cloud environments (see [the docs](http://snakemake.readthedocs.io/en/stable/executable.html) for details).

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.

### Step 4: Investigate results

After successful execution, you can create a self-contained report with all results via:

    snakemake --report report.html


Documentation
------------------
You can find the documentation [here](https://hoohm.github.io/dropSeqPipe/)

Future implementations
---------------------------
I'm actively seeking help to implement the points listed bellow. Don't hesitate to contact me if you wish to contribute.

* Create a sharing platform where quality plots/logs can be discussed and troubleshooted.
* Create a full html report for the whole pipeline
* Multiqc module for drop-seq-tools
* Implement an elegant "preview" mode where the pipeline would only run on a couple of millions of reads and allow you to have an approximated view before running all of the data. This would dramatically reduce the time needed to get an idea of what filters whould be used.
* 

I hope it can help you out in your single cell experiments!

Feel free to comment and point out potential improvements via [issues](https://github.com/Hoohm/dropSeqPipe/issues)


TODO
---------------------------------------------
* Add a mixed reference reference for testing purposes
* Finalize the parameters validation schema
* Make the debug feature a bit "cleaner". Deal with automatic naming of the debug variables
* Implement ddseq barcoding strategies
