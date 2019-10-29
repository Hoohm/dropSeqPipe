import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version

singularity:
    "shub://seb-mueller/singularity_dropSeqPipe:v04"

min_version("5.1.2")

#print(os.path.abspath(os.path.dirname(workflow.snakefile)))

# Load configuration files

try:
    configfile_path = config['configfile_path']
except:
    configfile_path = "config.yaml"
configfile: configfile_path


#Include the gtf biotypes yaml
configfile: config['META']['gtf_biotypes']

# Define a few variables to make them easier to reference
snakefile_root_path = os.path.abspath(os.path.dirname(workflow.snakefile))
ref_path = config['META']['reference-directory']
barcode_whitelist = config['FILTER']['barcode-whitelist']
results_dir = config['LOCAL']['results']
raw_data_dir = config['LOCAL']['raw_data']

# dropSeqPipe version
config['version'] = '0.5'
validate(config, schema=os.path.join(snakefile_root_path,"schemas","config.schema.yaml"))


# In order to deal with single species or mixed species experiment
# we define the same variables for each case.


#Define variables for mixed species experiments
if len(config['META']['species'].keys()) == 2:
    print('Running the pipeline for a mixed experiment')
    species_list = list(config['META']['species'])
    build_list = [
        config['META']['species'][species_list[0]]['build'],
        config['META']['species'][species_list[1]]['build']]
    release_list = [
        config['META']['species'][species_list[0]]['release'],
        config['META']['species'][species_list[1]]['release']]

    for species in config['META']['species']:
        release = '{}.{}'.format(
            config['META']['species'][species_list[0]]['release'],
            config['META']['species'][species_list[1]]['release'])
        build = '{}.{}'.format(
            config['META']['species'][species_list[0]]['build'],
            config['META']['species'][species_list[1]]['build'])
    species = 'mixed_{}_{}'.format(
        species_list[0],
        species_list[1])

#Define variables for single species experiments
elif len(config['META']['species'].keys()) == 1:
    species_list=list(config['META']['species'])
    species=species_list[0]
    release_list = [config['META']['species'][species]['release']]
    release=release_list[0]
    build_list = [config['META']['species'][species]['build']]
    build=build_list[0]
else:
    exit("Number of species in the config.yaml must be one or two. Exiting")

# Get sample names from samples.csv
samples = pd.read_table("samples.csv", sep=',').set_index("samples", drop=False)
validate(samples, schema=os.path.join(snakefile_root_path,"schemas","samples.schema.yaml"))
types=['read','umi']
# Get read_lengths from samples.csv
read_lengths = list(samples.loc[:,'read_length'])

wildcard_constraints:
    sample="({})".format("|".join(samples.index)),
    type="({})".format("|".join(types))


# Flexible ways to get the R1 and R2 files
def get_R1_files(wildcards):
    samples = [f for f in glob.glob("{}/*.fastq.gz".format(raw_data_dir)) if (re.search('R1', re.sub(wildcards.sample,'',f)) and re.search(wildcards.sample,f))]
    if len(samples)>1 & isinstance(samples,list):
        exit('Multiple read files for one sample. Please check file names or run snakemake -s rules/prepare.smk for multilane samples first.')
    if samples == []:
        exit('\tNo sample files found in the {}/ directory.\n\t\tPlease check that the path for the raw data is set properly in config.yaml'.format(raw_data_dir))
    return(samples)

def get_R2_files(wildcards):
    samples = [f for f in glob.glob("{}/*.fastq.gz".format(raw_data_dir)) if (re.search('R2', re.sub(wildcards.sample,'',f)) and re.search(wildcards.sample,f))]
    if len(samples)>1 & isinstance(samples,list):
        exit('Multiple read files for one sample. Please check file names or run snakemake -s rules/prepare.smk for multilane samples first.')
    if samples == []:
        exit('\tNo sample files found in the {} directory.\n\t\tPlease check that the path for the raw data is set properly in config.yaml'.format(raw_data_dir))
    return(samples)


if len(config['META']['species'].keys()) == 2:
    rule all:
        input:
            expand(
                ['{ref_path}/{species}_{build}_{release}/STAR_INDEX/SA_{read_length}/SA',
                #qc
                '{results_dir}/reports/fastqc_reads.html',
                '{results_dir}/reports/fastqc_barcodes.html',
                #fastqc_adapter
                'fastqc_adapter.tsv',
                #filter
                '{results_dir}/plots/adapter_content.pdf',
                '{results_dir}/reports/barcode_filtering.html',
                '{results_dir}/reports/RNA_filtering.html',
                '{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz',
                '{results_dir}/samples/{sample}/top_barcodes.csv',
                #mapping
                '{results_dir}/plots/knee_plots/{sample}_knee_plot.pdf',
                '{results_dir}/reports/star.html',
                '{results_dir}/plots/yield.pdf',
                '{results_dir}/samples/{sample}/Unmapped.out.mate1.gz',
                #splitting
                '{results_dir}/plots/barnyard/{sample}_genes.pdf',
                '{results_dir}/plots/barnyard/{sample}_transcripts.pdf'],
                    read_length=read_lengths,
                    sample=samples.index,
                    type=types,
                    results_dir=results_dir,
                    ref_path=config['META']['reference-directory'],
                    build=build,
                    release=release,
                    species=species),
            expand(
                ['{results_dir}/samples/{sample}/{species}/umi/matrix.mtx',
                '{results_dir}/samples/{sample}/{species}/read/matrix.mtx',
                '{results_dir}/plots/rna_metrics/{sample}_{species}_rna_metrics.pdf'],
                results_dir=results_dir,
                sample=samples.index,
                species=species_list)

elif len(config['META']['species'].keys()) == 1:
    rule all:
        input:
            #meta
            expand(
                ['{ref_path}/{species}_{build}_{release}/STAR_INDEX/SA_{read_length}/SA',
                #qc
                '{results_dir}/reports/fastqc_reads.html',
                '{results_dir}/reports/fastqc_barcodes.html',
                #filter
                '{results_dir}/plots/adapter_content.pdf',
                '{results_dir}/reports/barcode_filtering.html',
                '{results_dir}/reports/RNA_filtering.html',
                #mapping
                '{results_dir}/plots/knee_plots/{sample}_knee_plot.pdf',
                '{results_dir}/reports/star.html',
                '{results_dir}/plots/yield.pdf',
                '{results_dir}/samples/{sample}/Unmapped.out.mate1.gz',
                #extract
                '{results_dir}/plots/rna_metrics/{sample}_rna_metrics.pdf',
                '{results_dir}/summary/{type}/matrix.mtx',
                '{results_dir}/samples/{sample}/{type}/matrix.mtx',
                #merge
                '{results_dir}/plots/UMI_vs_counts.pdf',
                '{results_dir}/plots/UMI_vs_gene.pdf',
                '{results_dir}/plots/Count_vs_gene.pdf',
                '{results_dir}/summary/R_Seurat_objects.rdata',
                '{results_dir}/summary/barcode_stats_pre_filter.csv',
                '{results_dir}/summary/barcode_stats_post_filter.csv',
                '{results_dir}/plots/violinplots_comparison_UMI.pdf'],
                    read_length=read_lengths,
                    sample=samples.index,
                    type=types,
                    results_dir=results_dir,
                    ref_path=config['META']['reference-directory'],
                    build=build,
                    release=release,
                    species=species)
    rule download_meta:
        input:
            expand(
                ["{ref_path}/{species}_{build}_{release}/annotation.gtf",
                "{ref_path}/{species}_{build}_{release}/genome.fa"],
                    ref_path=config['META']['reference-directory'],
                    species=species_list,
                    release=release,
                    build=build)


rule qc:
    input:
        expand(
            ['{results_dir}/reports/fastqc_reads.html',
            '{results_dir}/reports/fastqc_barcodes.html',
            'fastqc_adapter.tsv'],
                results_dir=results_dir)

rule filter:
    input:
        expand(
            ['{results_dir}/plots/adapter_content.pdf',
            '{results_dir}/reports/barcode_filtering.html',
            '{results_dir}/reports/RNA_filtering.html',
            '{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz',
            '{results_dir}/samples/{sample}/top_barcodes.csv'],
                results_dir=results_dir,
                sample=samples.index)

rule map:
    input:
        expand(
            ['{results_dir}/plots/knee_plots/{sample}_knee_plot.pdf',
            '{results_dir}/reports/star.html',
            '{results_dir}/plots/yield.pdf',
            '{results_dir}/samples/{sample}/final.bam',
            '{results_dir}/samples/{sample}/Unmapped.out.mate1.gz'],
                sample=samples.index,
                results_dir=results_dir)

rule extract:
    input:
        expand(
            ['{results_dir}/plots/rna_metrics/{sample}_rna_metrics.pdf',
            '{results_dir}/summary/{type}/matrix.mtx',
            '{results_dir}/samples/{sample}/{type}/matrix.mtx.gz'],
                results_dir=results_dir,
                sample=samples.index,
                type=types)

rule split_species:
    input:
        expand(
            ['{results_dir}/samples/{sample}/{species}/barcodes.csv',
            '{results_dir}/plots/barnyard/{sample}_genes.pdf',
            '{results_dir}/plots/barnyard/{sample}_transcripts.pdf',
            '{results_dir}/samples/{sample}/{species}/unfiltered.bam'],
                sample=samples.index,
                species=config['META']['species'],
                results_dir=results_dir)


rule extract_species:
    input:
        expand(
            ['{results_dir}/samples/{sample}/{species}/{type}/matrix.mtx',
            '{results_dir}/plots/rna_metrics/{sample}_{species}_rna_metrics.pdf'],
                sample=samples.index,
                species=config['META']['species'],
                results_dir=results_dir,
                type=types)

rule merge:
    input:
        #merge
        expand(
            ['{results_dir}/plots/UMI_vs_counts.pdf',
            '{results_dir}/plots/UMI_vs_gene.pdf',
            '{results_dir}/plots/Count_vs_gene.pdf',
            '{results_dir}/summary/R_Seurat_objects.rdata',
            '{results_dir}/summary/barcode_stats_pre_filter.csv',
            '{results_dir}/summary/barcode_stats_post_filter.csv',
            '{results_dir}/plots/violinplots_comparison_UMI.pdf',
            '{results_dir}/summary/{type}/matrix.mtx'],
                results_dir=results_dir,
                type=types)

rule make_report:
    input:
        expand('{results_dir}/reports/publication_text.html', results_dir=results_dir)

if len(config['META']['species'].keys()) == 2:
    include: "rules/download_meta_mixed.smk"
if len(config['META']['species'].keys()) == 1:
    include: "rules/download_meta_single.smk"

include: "rules/generate_meta.smk"
include: "rules/fastqc.smk"
include: "rules/filter.smk"
include: "rules/cell_barcodes.smk"
include: "rules/map.smk"
include: "rules/extract_expression_single.smk"
include: "rules/split_species.smk"
include: "rules/extract_expression_species.smk"
include: "rules/merge.smk"
include: "rules/report.smk"
