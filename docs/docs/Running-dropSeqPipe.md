Context
-----------------------
The pipeline is supposed to be cloned once and then run on a folder containing your data as well as the configuration files. The run folder can contain multiple runs (aka batches) as you can easily add new samples when recieving new data and run the same commands. This will simply run the pipeline on the newly added data and recreate reports as well as plots containing all the samples.

Example: You run 2 biological conditions with 2 replicates. This makes up for 4 samples. Assume a simple dropseq protocol with only human cells.
1. You sequence the data and recieve the 8 files (two files per sample) and set up the pipeline
2. You run the pipeline with the command: `snakemake --use-conda --cores N --directory WORKING_DIR` `N` being the number of cores available and `WORKING_DIR` being the folder containing your data.
3. You see that there is an issue with the protocol and you modify it
4. You create a new set of libraries and sequence them (same 2x2 design)
5. You add the new files in the data folder of `WORKING_DIR` and edit the samples.csv to add missing samples.
6. You run the pipeline as you did the first time `snakemake --use-conda --cores N --directory WORKING_DIR`
7. This will run the new samples only and recreate the reports as well as the BC_drop and yield plots.
8. It is now easy to compare the impact of your change in the procotol

Working dir folder preparation
----------------
The raw data from the sequencer should be stored in the `data` folder of `WORKING_DIR` folder like this:
```
/path/to/your/WORKING_DIR/
| -- data/
| -- -- sample1_R1.fastq.gz
| -- -- sample1_R2.fastq.gz
| samples.csv
| config.yaml
| barcodes.csv
```
*Note: In DropSeq or ScrbSeq you expect a paired sequencing. R1 will hold the information of your barcode and UMI, R2 will hold the 3' end of the captured mRNA.*


Once everything is in place, you can run the pipeline using the normal snakemake commands.

Running the pipeline (TLDR version)
----------------------------

For a simple single cell run you only need to run: `snakemake --cores N --use-conda --directory WORKING_DIR`
This will run the whole pipeline and use the X number of cores you gave to it.


Running the pipeline
---------------------------------

I highly recommend to take a [look at the options](http://snakemake.readthedocs.io/en/latest/) that are available since I won't cover everything here.


Modes
------------------------------
You have two main ways to run the pipeline.

You can either just run `snakemake --use-conda --directory WORKING_DIR` in the root folder containing your experiment and it will run everything without stopping.

You can also run each step separately. The main advantage of the second way is that you are able to fine tune your parameters based on the results of fastqc, filtering, mapping quality, etc...
I would suggest using the second approach when you work on a new protocol and the first one when you are confident of your parameters.

There are seven different modes available and to run one specifically you need to call the mode.

Example: To run the `qc` mode: `snakemake --cores 8 qc --use-conda --directory WORKING_DIR`
You can also run multiple modes at the same time if you want: `snakemake --cores 8 qc filter --use-conda --directory WORKING_DIR`

### Single species:
* `meta`: Generates all the subsequent references files and STAR index needed to run the pipeline. You can run this alone if you just want to create the meta-data file before running a new set of data.
* `qc`: Creates fastqc reports of your data.    
* `filter`: Go from sample_R1.fastq.gz to sample_filtered.fastq.gz ready to be mapped to the genome.  This step filters out your data for low quality reads and trims adapter you provided in the FILTER section.
* `map`: Go from sample_filtered.fastq.gz to the sample_final.bam read to extract the expression data. This maps the data to the genom.
* `extract`: Extract the expression data. You'll get a umi and a count expression matrix from your whole experiment.

### Mixed species
The next two step are only necessary for mixed experiment. Note that you can't just run `snakemake --directory WORKING_DIR --use-conda` for a mixed experiment but you have to run `snakemake extract-species --directory WORKING_DIR` instead.
* `split_species`: Poduces the species plot and the species specific barcode list for extract-species.
* `extract_species`: Extract the expression data. You'll get a umi and a count expression matrix from your whole experiment/specie in the summary folder.

If you want all the plots for mixed experiment, you have to call all the modes:  `snakemake --cores 8 qc filter map split_species extract_species --use-conda --directory WORKING_DIR`

Barcode whitelist
---------------------
In protocols such as SCRBseq, the expected barcodes sequences are known. This pipeline also does allow the use of known barcodes instread of a number of expected cells.
In order to use this functionnality you just need to add a `barcodes.csv` file at the root of the working dir and run the pipeline as usual.

Advanced options
-------------------
If you have some specific adapters that are not present by default in the ones in the `templates` folder, you can add whatever adapters you want to trim (as many as you need) following the fasta syntax in any fasta files proposed. Be sure to specify it properly in the config.yaml at `FILTER:IlluminaClip`

Further options
---------------------
* `--cores N` Use this argument to use X amunt of cores available.
* `--notemp` Use this to not delete all the temporary files. Without this option, only files between steps are kept. Use this option if you are troobleshooting the pipeline or you want to analyze in between files yourself.
* `--dryrun` Use this to check out what is going to run if you run your command. This is nice to check for potential missing files.


Running on clusters
----------------------------------
There is a file in the `templates` called `cluster.yaml`. This can be used to modify ressources needed for your data. I generally recommand moving the file to the root of the folder so that it doesn't get replaced by updates.

Bellow is an example of running on a cluster using the template file `cluster.yaml` on SLURM.

```
snakemake --cluster 'sbatch -n {cluster.n}  -t {cluster.time} --clusters=CLUSTERNAME --output={cluster.output}' --jobs N --cluster-config cluster.yaml --use-conda --local-cores C
```

* N: is the number of jobs you are allowed to run at the same time
* C: is the local-cores of the host machine. A few simple rules are gonna be run locally (not sent to nodes) because they are not that heavy (mostly plotting)
* CLUSTERNAME: the name of the cluster you want to use


Folder Structure
-----------------------
This is the folder structure you get in the end:
```
/path/to/your/WORKING_DIR/
| -- data/
| -- logs/
| -- -- cluster/
| -- plots/
| -- reports/
| -- summary/
| Snakefile
| samples.csv
| config.yaml
| barcodes.csv
| .snakemake/
```

* `data/` Contains all your samples as well as the intermediary files
* `logs/` Contains all the logfiles generated by the pipeline
* `logs/cluster` Contains all the logfiles generated by the cluster
* `plots/` Contains all the plots generated by the pipeline
* `reports/` Contains all the reports generated by the pipeline
* `summary/` Contains all the files you might use for downstream analysis (contains barcodes selected per sample per species, final umi/counts expression matrix)
* `samples.csv` File containing sample details
* `config.yaml` File containing pipeline parameters as well as system parameters
* `.snakemake/` Folder that contains all the environements created for the run as well as a lot of other things that I don't know about.
