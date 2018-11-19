import pandas as pd
import os
import re
import glob

ruleorder: extend_barcode_whitelist > extend_barcode_top

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
    samples = [f for f in glob.glob("data/*.fastq.gz") if (re.search('R1', re.sub(wildcards.sample,'',f)) and re.search(wildcards.sample,f))]
    if len(samples)!=1:
        exit('Multiple read files for one sample. Please check file names or run snakemake -s rules/prepare.smk for multilane samples first.')
    return(samples)

def get_R2_files(wildcards):
    samples = [f for f in glob.glob("data/*.fastq.gz") if (re.search('R2', re.sub(wildcards.sample,'',f)) and re.search(wildcards.sample,f))]
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
        #expand('data/{sample}/salmon/mapping.tsv', sample=samples.index),
        #extract
        'summary/umi/expression.mtx',
        'summary/reads/expression.mtx',
        'summary/R_Seurat_objects.rdata',
        expand('plots/rna_metrics/{sample}_rna_metrics.pdf', sample=samples.index),
        #expand('data/{sample}/{type}/expression_matrix.mtx', sample=samples.index, type=types),

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
        'plots/yield.pdf',        
        'reports/star.html'
        
rule extract:
    input:
        expand('plots/rna_metrics/{sample}_rna_metrics.pdf', sample=samples.index),
        expand('plots/{sample}_{type}_expression.long', sample=samples.index, type=types),
        

rule split_species:
    input:
        expand('plots/knee_plots/{sample}_knee_plot.pdf', sample=samples.index),
        expand('data/{sample}/{species}/barcodes.csv', sample=samples.index, species=config['META']['species']),
        expand('plots/barnyard/{sample}_genes.pdf', sample=samples.index),
        expand('plots/barnyard/{sample}_transcripts.pdf', sample=samples.index),
        expand('data/{sample}/{species}/unfiltered.bam', sample=samples.index, species=config['META']['species'])


rule extract_species:
    input:
        expand('data/{sample}/{species}/umi_expression_matrix.txt', sample=samples.index, species=config['META']['species']),
        expand('data/{sample}/{species}/counts_expression_matrix.txt', sample=samples.index, species=config['META']['species']),
        expand('summary/Experiment_{species}_counts_expression_matrix.tsv', species=config['META']['species']),
        expand('summary/Experiment_{species}_umi_expression_matrix.tsv', species=config['META']['species']),
        expand('plots/rna_metrics/{sample}_{species}_rna_metrics.pdf', sample=samples.index, species=config['META']['species'])
        
rule merge:
    input:
        #merge
        expand('summary/{sample}/{type}/expression.mtx', sample=samples.index, type=types),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv',
        'plots/violinplots_comparison_UMI.pdf',
        'plots/UMI_vs_counts.pdf',
        'plots/UMI_vs_gene.pdf',
        'plots/Count_vs_gene.pdf',
        'summary/R_Seurat_objects.rdata'



include: "rules/generate_meta.smk"
include: "rules/fastqc.smk"
include: "rules/filter.smk"
include: "rules/cell_barcodes.smk"
include: "rules/map.smk"

include: "rules/extract_expression_single.smk"

if(config['META']['mixed']):
    include: "rules/split_species.smk"
    include: "rules/extract_expression_species.smk"

include: "rules/merge.smk"
