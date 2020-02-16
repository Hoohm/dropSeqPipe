Welcome
------------------------------

Welcome to the documentation of dropSeqPipe v`0.5`.

DropSeqPipe (dSP)is a pipeline for the analysis of single cell data. 
The pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It allows to go from raw data of a single cell experiment to the final count matrix including extensive QC readouts. 
dSP provides an easy framework to reproduce and compare different experiments with different parameters using STAR as an aligner. It is usable for any single cell protocol using paired-end reads where the first read holds the cell and UMI barcodes and the second read the transcript. 
Below is a non-exhaustive list of compatible protocols and brands:

* Drop-Seq
* SCRB-Seq
* 10x Genomics
* DroNc-seq
* Dolomite Bio ([Nadia Instrument](https://www.dolomite-bio.com/product/nadia-instrument/))
    
This package is trying to be as user friendly as possible. Even though it requires command line executions, the hope is that non-bioinformatician would be able to use it. 

Note: A short demonstration on how to to install and use dSP can be found in [this webinar by Sebastian Mueller](https://www.youtube.com/watch?v=4bt-azBO-18).
