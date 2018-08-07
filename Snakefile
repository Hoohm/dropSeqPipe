import pandas as pd
import os
import re

# Load configuration file
configfile: "config.yaml"

# Get sample names from samples.csv
samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)

# Get read_lengths from samples.csv
read_lengths = list(samples.loc[:,'read_length'])

# Constraint sample names wildcards
wildcard_constraints:
    sample="({})".format("|".join(samples.index))

# Create reference files prefixes
reference_prefix = os.path.join(config['META']['reference-directory'], re.split(".fasta|.fa",config['META']['reference-file'])[0])
annotation_prefix = os.path.join(config['META']['reference-directory'],config['META']['annotation-file'].split('.gtf')[0])
reference_file = os.path.join(config['META']['reference-directory'], config['META']['reference-file'])
annotation_file = os.path.join(config['META']['reference-directory'], config['META']['annotation-file'])
annotation_reduced_file = os.path.join(config['META']['reference-directory'],'.'.join([config['META']['annotation-file'].split('.gtf')[0],'reduced','gtf']))
star_index_prefix = os.path.join(config['META']['reference-directory'],'STAR_INDEX/SA')

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
        'reports/filter.html',
        #mapping
        expand('plots/{sample}_knee_plot.pdf', sample=samples.index),
        'reports/star.html',
        'plots/yield.pdf',
        'plots/UMI_vs_counts.pdf',
        'plots/UMI_vs_gene.pdf',
        'plots/Count_vs_gene.pdf',
        'summary/R_Seurat_objects.rdata',
        #extract
        expand('logs/dropseq_tools/{sample}_umi_per_gene.tsv', sample=samples.index),
        expand('plots/{sample}_rna_metrics.pdf', sample=samples.index),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv'


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
        expand('data/{sample}/filtered.fastq.gz', sample=samples.index),
        'reports/filter.html',
        'plots/adapter_content.pdf'
        
rule map:
    input:  
        expand('data/{sample}/final.bam', sample=samples.index),
        expand('logs/dropseq_tools/{sample}_hist_out_cell.txt', sample=samples.index),
        expand('plots/{sample}_knee_plot.pdf', sample=samples.index),
        'reports/star.html',
        'plots/violinplots_comparison_UMI.pdf',
        # 'plots/UMI_vs_counts.html',
        'plots/UMI_vs_counts.pdf',
        # 'plots/UMI_vs_gene.html',
        'plots/UMI_vs_gene.pdf',
        # 'plots/Count_vs_gene.html',
        'plots/Count_vs_gene.pdf',
        'summary/R_Seurat_objects.rdata',
        'plots/yield.pdf'
        
rule extract:
    input:
        expand('logs/dropseq_tools/{sample}_umi_per_gene.tsv', sample=samples.index),
        expand('plots/{sample}_rna_metrics.pdf', sample=samples.index),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv'
        

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
        expand('plots/{species}/{sample}_rna_metrics.pdf', sample=samples.index, species=config['META']['species'])
        


include: "rules/generate_meta.smk"
include: "rules/fastqc.smk"
include: "rules/filter.smk"
include: "rules/map.smk"
include: "rules/extract_expression_single.smk"
include: "rules/split_species.smk"
include: "rules/extract_expression_species.smk"
