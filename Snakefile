import pandas as pd
import os
import re
import glob

# Load configuration file
configfile: "config.yaml"

# Get sample names from samples.csv
samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)
types=['counts','umi']
# Get read_lengths from samples.csv
read_lengths = list(samples.loc[:,'read_length'])

# Constraint sample names wildcards
wildcard_constraints:
    sample="({})".format("|".join(samples.index))


def get_R1_files(wildcards):
    samples = [f for f in glob.glob("data/*.fastq.gz") if (re.search('R1', f) and re.search(wildcards.sample,f))]
    if len(samples)!=1:
        exit('Multiple read files for one sample. Please check file names or run snakemake -s rules/prepare.smk for multilane samples first.')
    return(samples)

def get_R2_files(wildcards):
    samples = [f for f in glob.glob("data/*.fastq.gz") if (re.search('R2',f) and re.search(wildcards.sample,f))]
    if len(samples)!=1:
        exit('Multiple read files for one sample. Please check file names or run snakemake -s rules/prepare.smk for multilane samples first.')
    return(samples)



# Create reference files prefixes
reference_prefix = os.path.join(config['META']['reference-directory'], re.split(".fasta|.fa",config['META']['reference-file'])[0])
annotation_prefix = os.path.join(config['META']['reference-directory'],config['META']['annotation-file'].split('.gtf')[0])
reference_file = os.path.join(config['META']['reference-directory'], config['META']['reference-file'])
annotation_file = os.path.join(config['META']['reference-directory'], config['META']['annotation-file'])
annotation_reduced_file = os.path.join(config['META']['reference-directory'],'.'.join([config['META']['annotation-file'].split('.gtf')[0],'reduced','gtf']))
star_index_prefix = os.path.join(config['META']['reference-directory'],'STAR_INDEX/SA')
salmon_index = os.path.join(config['META']['reference-directory'], 'salmon','index')
# Get barcode length
starttrim_length = config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1




rule all:
    input:
        #meta
        '{}.refFlat'.format(annotation_prefix),
        '{}.reduced.gtf'.format(annotation_prefix),
        '{}.dict'.format(reference_prefix),
        '{}.rRNA.intervals'.format(reference_prefix),
        expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths),
        #qc
        'reports/fastqc_reads.html',
        'reports/fastqc_barcodes.html',
        #filter
        'plots/adapter_content.pdf',
        'reports/barcode_filtering.html',
        'reports/RNA_filtering.html',
        #mapping
        expand('plots/knee_plots/{sample}_knee_plot.pdf', sample=samples.index),
        'reports/star.html',
        'plots/yield.pdf',
        'plots/UMI_vs_counts.pdf',
        'plots/UMI_vs_gene.pdf',
        'plots/Count_vs_gene.pdf',
        'summary/R_Seurat_objects.rdata',
        #expand('data/{sample}/salmon/mapping.tsv', sample=samples.index),
        #extract
        expand('logs/dropseq_tools/{sample}_umi_per_gene.tsv', sample=samples.index),
        expand('plots/rna_metrics/{sample}_rna_metrics.pdf', sample=samples.index),
        #expand('data/{sample}/{type}/expression_matrix.mtx', sample=samples.index, type=types),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv'



rule test:
    input:
        expand('data/{sample}_R1.fastq.gz', sample=samples.index),
        expand('data/{sample}_R2.fastq.gz', sample=samples.index)


rule meta:
    input:
        '{}.refFlat'.format(annotation_prefix),
        '{}.reduced.gtf'.format(annotation_prefix),
        '{}.dict'.format(reference_prefix),
        '{}.rRNA.intervals'.format(reference_prefix),
        expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths)

rule qc:
    input:
        'reports/fastqc_reads.html',
        'reports/fastqc_barcodes.html'

rule filter:
    input:
        expand('data/{sample}/trimmmed_repaired_R1.fastq.gz', sample=samples.index),
        'reports/barcode_filtering.html',
        'reports/RNA_filtering.html',
        'plots/adapter_content.pdf'
        
rule map:
    input:  
        expand('data/{sample}/final.bam', sample=samples.index),
        expand('logs/dropseq_tools/{sample}_hist_out_cell.txt', sample=samples.index),
        expand('plots/knee_plots/{sample}_knee_plot.pdf', sample=samples.index),
        #expand('data/{sample}/salmon/mapping.tsv', sample=samples.index),
        'reports/star.html'
        
rule extract:
    input:
        expand('logs/dropseq_tools/{sample}_umi_per_gene.tsv', sample=samples.index),
        expand('plots/rna_metrics/{sample}_rna_metrics.pdf', sample=samples.index),
        #expand('data/{sample}/{type}/expression_matrix.mtx', sample=samples.index, type=types),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv',
        'plots/violinplots_comparison_UMI.pdf',
        'plots/UMI_vs_counts.pdf',
        'plots/UMI_vs_gene.pdf',
        'plots/Count_vs_gene.pdf',
        'summary/R_Seurat_objects.rdata',
        'plots/yield.pdf'
        

rule split_species:
    input:
        expand('summary/{species}/{sample}_barcodes.csv', sample=samples.index, species=config['META']['species']),
        expand('plots/{sample}_species_plot_genes.pdf', sample=samples.index),
        expand('plots/{sample}_species_plot_transcripts.pdf', sample=samples.index),
        expand('data/{species}/{sample}_unfiltered.bam', sample=samples.index, species=config['META']['species'])


rule extract_species:
    input:
        expand('summary/{species}/{sample}_umi_expression_matrix.txt', sample=samples.index, species=config['META']['species']),
        expand('summary/{species}/{sample}_counts_expression_matrix.txt', sample=samples.index, species=config['META']['species']),
        expand('logs/{species}/{sample}_umi_per_gene.tsv', sample=samples.index, species=config['META']['species']),
        expand('summary/Experiment_{species}_counts_expression_matrix.tsv', species=config['META']['species']),
        expand('summary/Experiment_{species}_umi_expression_matrix.tsv', species=config['META']['species']),
        expand('plots/{species}/rna_metrics/{sample}_rna_metrics.pdf', sample=samples.index, species=config['META']['species'])
        


include: "rules/generate_meta.smk"
include: "rules/fastqc.smk"
include: "rules/filter.smk"
include: "rules/map.smk"
include: "rules/extract_expression_single.smk"
include: "rules/split_species.smk"
include: "rules/extract_expression_species.smk"