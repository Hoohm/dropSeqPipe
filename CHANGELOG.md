# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).


## [0.23]
### Changed
- pre_align steps will output a fastq.gz instead of a fastq file.
- `fastqc.R` is now compatible with paired and single end data.
- Changed a few options in `GLOBAL` for `UMI` and `Cell_barcodes` options. Now possible to change filtering settings. See README
- STAR logs have been stripped of the `STAR` string. This is to allow for better compatibility with [multiqc](https://github.com/ewels/MultiQC/)
- Removed `fastqc` folder and moved items to `logs` folder. Grouping all logs files for better [multiqc](https://github.com/ewels/MultiQC/) compatibility.
- Changed `generate_meta` to `generate-meta` for keeping similar syntax between modes.
- Added seperate log files for stats and summary in the DetectBeadSynthesisErrors.

### Added
- Bulk data single or paired end data processing capacity.
- Started a wiki with a FAQ


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