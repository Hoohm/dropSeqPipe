import pandas as pd
import os

configfile: "config.yaml"
samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)

wildcard_constraints:
    sample="({})".format("|".join(samples.index))


reference_prefix = os.path.join(config['META']['reference-directory'], config['META']['reference-file'].split('.fasta')[0])
annotation_prefix = os.path.join(config['META']['reference-directory'],config['META']['annotation-file'].split('.gtf')[0])
reference_file = os.path.join(config['META']['reference-directory'], config['META']['reference-file'])
annotation_file = os.path.join(config['META']['reference-directory'], config['META']['annotation-file'])
annotation_reduced_file = os.path.join(config['META']['reference-directory'],'.'.join([config['META']['annotation-file'].split('.gtf')[0],'reduced','gtf']))
star_index_prefix = os.path.join(config['META']['reference-directory'],'STAR_INDEX/SA')


starttrim_length = config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1

read_lengths = list(samples.loc[:,'read_length'])


rule all:
    input:
        #meta
        '{}.refFlat'.format(annotation_prefix),
        '{}.reduced.gtf'.format(annotation_prefix),
        '{}.dict'.format(reference_prefix),
        '{}.rRNA.intervals'.format(reference_prefix),
        expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths),
        #qc
        expand('logs/{sample}_R1_fastqc.html', sample=samples.index),
        expand('logs/{sample}_R2_fastqc.html', sample=samples.index),
        'reports/fastqc_reads.html',
        'reports/fastqc_barcodes.html',
        'reports/fastqc_reads_data/multiqc_general_stats.txt',
        #filter
        expand('data/{sample}_filtered.fastq.gz', sample=samples.index),
        expand('plots/{sample}_polya_trimmed.pdf', sample=samples.index),
        expand('plots/{sample}_start_trim.pdf', sample=samples.index),
        expand('plots/{sample}_CELL_dropped.pdf', sample=samples.index),
        expand('plots/{sample}_UMI_dropped.pdf', sample=samples.index),
        'plots/BC_drop.pdf',
        'reports/filter.html',
        #mapping
        expand('data/{sample}_final.bam', sample=samples.index),
        expand('logs/{sample}_hist_out_cell.txt', sample=samples.index),
        expand('plots/{sample}_knee_plot.pdf', sample=samples.index),
        'reports/star.html',
        'plots/yield.pdf',
        #extract
        expand('logs/{sample}_umi_per_gene.tsv', sample=samples.index),
        expand('plots/{sample}_rna_metrics.pdf', sample=samples.index),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv',
        'reports/combined.html'
        
rule meta:
    input:
        '{}.refFlat'.format(annotation_prefix),
        '{}.reduced.gtf'.format(annotation_prefix),
        '{}.dict'.format(reference_prefix),
        '{}.rRNA.intervals'.format(reference_prefix),
        expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths)

rule qc:
    input:
        expand('logs/{sample}_R1_fastqc.html', sample=samples.index),
        expand('logs/{sample}_R2_fastqc.html', sample=samples.index),
        'reports/fastqc_reads.html',
        'reports/fastqc_barcodes.html',
        'reports/fastqc_reads_data/multiqc_general_stats.txt'

rule filter:
    input:
        expand('data/{sample}_filtered.fastq.gz', sample=samples.index),
        expand('plots/{sample}_polya_trimmed.pdf', sample=samples.index),
        expand('plots/{sample}_start_trim.pdf', sample=samples.index),
        expand('plots/{sample}_CELL_dropped.pdf', sample=samples.index),
        expand('plots/{sample}_UMI_dropped.pdf', sample=samples.index),
        'reports/filter.html',
        'plots/BC_drop.pdf'
        
rule map:
    input:  
        expand('data/{sample}_final.bam', sample=samples.index),
        expand('logs/{sample}_hist_out_cell.txt', sample=samples.index),
        expand('plots/{sample}_knee_plot.pdf', sample=samples.index),
        'reports/star.html',
        'plots/yield.pdf'
        
rule extract:
    input:
        expand('logs/{sample}_umi_per_gene.tsv', sample=samples.index),
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
        

rule multiqc_all:
    input: 
        fastqc = expand('logs/{sample}_R1_fastqc.html',sample=samples.index),
        trimmomatic=expand('logs/{sample}_trimlog.txt',sample=samples.index),
        star=expand('data/{sample}/Log.final.out',sample=samples.index)
    output: 'reports/combined.html'
    shell:
        """multiqc logs/ -m fastqc -m trimmomatic -m star -o reports/ -n combined.html -f"""

include: "rules/generate_meta.smk"
include: "rules/fastqc.smk"
include: "rules/filter.smk"
include: "rules/map.smk"
include: "rules/extract_expression_single.smk"
include: "rules/split_species.smk"
include: "rules/extract_expression_species.smk"