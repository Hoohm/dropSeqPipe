# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).


## [0.5]
### Added
- Singularity usage. Try out the `--use-singularity` option instead of `--use-conda`

### Changed 
- Lots off small bugfixes


## [0.4.1]
### Added
- samples.csv and config.yaml schema validation. This will help users fix missing values.
- DetectBeadSubstitutionErrors was added in the mapping steps.

### Changed
- Minimum read length after trimming is now the index of the end of the UMI
- dropSeqPipe can now run with a docker image if you use the `--use-singularity` option. This should help people with package issues and different linux setups. You need to have installed singularity system wide to use this option.


## [0.4] - 2018-12-19
### Added
- Top barcode detection using [umi-tools](https://github.com/CGATOxford/UMI-tools) based on number of expected cells.
- Genome reference and annotation automatically downloaded now base on build and release number from configuration file.
- On the fly detection of mixed experiment.
- **beta**: Generation of a report for publication describing tools used in each steps. run `make_report` after the preprocessing is done to get `reports/publication_text.html`. This is a really early stage. Feel free to suggest PR for text modifications.
- Raw data, results, reference are now independent from the working dir and can be chosen via the configuration file.
- dropseq_tools v2.0 implemented. This opens up new options such as choosing which locus to use for gene counting. See configuration file.
- Possibility to edit which biotypes are selected from the annotation file via a gtf_biotypes.yaml file provided.
- Cell barcodes are now corrected. One hamming distance for known/given whitelists, graphbased correction based on umi-tools for unknown lists. Those corrections are written in the bam files. This makes final bam files compatible for other tools using the XC/XM bam TAGS.
- UMI are now also corrected based on dropseq_tools v2.0.
- Possibility to choose SENSE, ANTISENSE or BOTH for read counting.
- Adapter content for R1 and R2 have now their own plot, `adapter_content.pdf`.
- New plot called `yield.pdf` makes a summary of total reads and how they are distributed among filtered, trimmed, mapped, etc.
- Configuration file has now a CONTACT section providing a field for a person and a contact e-mail address.

### Changed
- Expression matrices output are now sparse (mtx format). This will decrease the size of the output and loading time for downstream analysis.
- Logfiles, plots and samples output are now grouped together in folders by category. This should make browsing results easier.
- Fixed most of the packages versions.
- Summary plots and Seurat object are now in the `all` rule and will be created by default.

### Removed
- Merging of species expression accross samples. Since the mixed experiments are mostly used to test out the doublet rate of a platform and not for downstream analysis, this last part has not been updated. Single expression matrices are still there.
- Cell barcodes dropped, umi barcodes dropped, starttrim and polyA trim plots are now gone. BC_drop is also removed. Replacements are adapter_content and yield plots.
- Quality trimming via dropseq_tools has been removed and is now down by cutadapt. Those modifications decrease the running time of the pipeline.


## [0.32]
### Added
- Documentation generated from the markdown files directly on travis-ci.


## [0.31a]
### Changed
- fix on species plot.
- fix on rule STAR_align adding now unmapped read to a fastq file.

### Added
- Added travis integration. The pipeline is now automatically getting tested when updated and when pull requests are proposed.
- There is now a small git submodule in .test which will provide a sampled file for testing the pipeline on travis-ci.

### Removed
- `environment.yaml` has been removed. Youjust have to install snakemake now instead of activating the env.

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


## [0.3]
### Changed
- Complete overhaul of how the pipeline is organized to follow the structure proposed for snakemake-workflows. This will allow ease of deployement on any platform having conda installed. It will also help to run on clusters.
- The way to call the pipeline is now simplified. Changes are shown in the [WIKI](https://github.com/Hoohm/dropSeqPipe/wiki/)
- Dependency to Drop-seq-tools updated from version 1.12 to 1.13
- Full compatibility with barcode whitelist. Makes it easier to use for SCRBseq protocols or whitelist from other source (UMI-tools).
- Modified cell and UMI drop plots in order to reflect the option chosen. See [plots](https://github.com/Hoohm/dropSeqPipe/wiki/Plots)

### Removed
- Bulk sequencing compatiblity.
- Fastqc and STAR logs plots are removed and replaced by multiqc.
- Automatic determination of STAMPS via knee_plot. Please use an estimated number of cells as the main threshold and filter in downstream analysis for other parameters such as high number of mitochondrial genes.
- `MinCellFraction` entry in config.yaml. This parameter wasn't adding much value and was confusing.
- Base frequency plot has been removed. This will come back with autodetermination of the STAMPS.

### Added
- Wrapper for Drop-seq tools. Makes it easier to switch temp folder and choose maximum memory heap.
- More parameters for STAR exposed. See [WIKI](https://github.com/Hoohm/dropSeqPipe/wiki/)

## [0.24]
### Changed
- All the QCplots are now generated inside the snakefiles. No more `generate-plots` mode.


## [0.23a]
### Changed
- Will now allow you to run `generate-meta` without having a `config.yaml` file in the reference foder.
- Changed the code for Cell and UMI barcode quality drop (per sample and overall). There was an error in the code not givint the right amount of dropped reads. Updated the images on the wiki accordingly.
- Fixed the setup where r2py was called before getting installed.
- Big change in the mapping. From now on the STAR index will be done without a GTF file. This allows to change the overhang option on the fly for each sample based on the mean read length. This also opens up 2-pass mapping. You will have to regenerate your index for it to work.
- Changed `generate_meta` in order to fit the new STAR index without a GTF. You now have to give the path to the GTF file in the config.yaml

### Added
- `min_count_per_umi` in the `config.yaml` to decide how many times a Gene - UMI has to be found to be counted as one.


## [0.23]
### Changed
- pre_align steps will output a fastq.gz instead of a fastq file.
- `fastqc.R` is now compatible with paired and single end data.
- Changed a few options in `GLOBAL` for `UMI` and `Cell_barcodes` options. Now possible to change filtering settings. See [WIKI](https://github.com/Hoohm/dropSeqPipe/wiki/Create-config-files)
- STAR logs have been stripped of the `STAR` string. This is to allow for better compatibility with [multiqc](https://github.com/ewels/MultiQC/)
- Removed `fastqc` folder and moved items to `logs` folder. Grouping all logs files for better [multiqc](https://github.com/ewels/MultiQC/) compatibility.
- Changed `generate_meta` to `generate-meta` for keeping similar syntax between modes.
- Added seperate log files for stats and summary in the DetectBeadSynthesisErrors.
- Moved part of the `README`to the wiki.
- Changed the name of the first expression matrix extracted before the species plot to `unfiltered_expression.`


### Added
- You can now run Bulk Single or paired end RNAseq data.
- Started a wiki with a FAQ
- Added options in `GLOBAL` config.yaml. You can now choose a range of options for UMI and Barcode filtering. please refer to the wiki for more information.
- Support for [MultiQC](https://github.com/ewels/MultiQC/). MultiQC is a great way of summarising all of the logs from your experiment. As of today it supports 46 different modules (such as fastqc, trimmomatic, STAR, etc...) The `generate-plots` mode now produces a `multiqc_report.html` file in the plots folder.
- New plot! BCDrop.pdf is a new plot showing you how many barcode and UMIs you dropped from the raw data before aligning. This helps to track how many samples you might loose because of low quality reads in the barcoding.

## [0.22]
### Changed
- all `subprocess.call` replaced by `shell` from snakemake
- STAR aligner now not limited to 8 cores or threads but will use the maximum number provided in the local.yaml file
- Name from dropSeqPip to dropSeqPipe
- Fixed a bug where all stage1 steps used the same summary file. Now BC tagging, UMI tagging, starting trim and polyA trim have different summary files
- extract-expression now merges all the samples final count matrix into one per run (folder)
- Fixed a bug where the amount of total reads on the knee-plot was overinflated.
- Changed `knee-plot` mode to `generate-plots`.

### Added
- Temp files have been added in the pipeline. You can turn this off by using the `--notemp` option
- fastqc mode now available. Generates fastqc reports plus summary plots
- Summary file and plot for fastqc and STAR logs
- Missing R packages should install automatically now. No need to install them beforehand. Report any problem plz
- `GLOBAL` values in the config files are now available. They allow to change UMI and BC ranges as well as mismatches for STAR aligner
- Added a new mode: generate_meta. This allows to create all the metadata files needed for the pipeline. You just need a folder with a genome.fa and an annotation.gtf

## [0.21]
### Added
- Changelog file to track changes
- --rerun option to force a rerun
- Multiple steps allowd now

## [0.2] - 2017-03-14
### Changed
- The pipeline is now a python package being called as an executable
- Went from json to yaml for config files

### Added
- setup.py and dependencies
- Species plot available

### Removed
- primer handling, went to default: AAGCAGTGGTATCAACGCAGAGTAC


## [0.1] - 2017-02-13
### First release
- Allows for preprocessing, alignement with STAR, post align processing until knee-plot