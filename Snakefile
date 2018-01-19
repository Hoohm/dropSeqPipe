import pandas as pd
import os

configfile: "config.yaml"

reference_prefix = os.path.join(config['META']['reference_folder'], config['META']['reference_file'].split('.fasta')[0])
annotation_prefix = os.path.join(config['META']['reference_folder'],config['META']['annotation_file'].split('.gtf')[0])
reference_file = os.path.join(config['META']['reference_folder'], config['META']['reference_file'])
annotation_file = os.path.join(config['META']['reference_folder'], config['META']['annotation_file'])
annotation_reduced_file = os.path.join(config['META']['reference_folder'],'.'.join([config['META']['annotation_file'].split('.gtf')[0],'reduced','gtf']))
star_index_prefix = os.path.join(config['META']['reference_folder'],'STAR_INDEX/SA')


starttrim_length = config['FILTER']['Cell_barcode']['end'] - config['FILTER']['Cell_barcode']['start'] + 1

samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)
read_lengths = list(samples.loc[:,'read_length'])


rule all:
    input:
        #meta
        expand('{annotation_prefix}.refFlat',annotation_prefix=annotation_prefix),
        expand('{annotation_prefix}.reduced.gtf',annotation_prefix=annotation_prefix),
        expand('{reference_prefix}.dict',reference_prefix=reference_prefix),
        expand('{reference_prefix}.rRNA.intervals',reference_prefix=reference_prefix),
        expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths),
        #qc
        expand('logs/{sample}_R1_fastqc.html', sample=samples.index),
        'reports/fastqc.html',
        #filter
        expand('data/{sample}_filtered.fastq.gz', sample=samples.index),
        expand('plots/{sample}_polya_trimmed.pdf', sample=samples.index),
        expand('plots/{sample}_start_trim.pdf', sample=samples.index),
        expand('plots/{sample}_CELL_dropped.pdf', sample=samples.index),
        expand('plots/{sample}_UMI_dropped.pdf', sample=samples.index),
        'reports/filter.html',
        'plots/BC_drop.pdf',
        'plots/png/BC_drop.png',
        #mapping
        expand('data/{sample}_final.bam', sample=samples.index),
        'reports/star.html',
        expand('logs/{sample}_hist_out_cell.txt', sample=samples.index),
        #extract
        expand('logs/{sample}_umi_per_gene.tsv', sample=samples.index),
        expand('plots/{sample}_rna_metrics.pdf', sample=samples.index),
        expand('plots/png/{sample}_rna_metrics.png', sample=samples.index),
        expand('plots/{sample}_knee_plot.pdf', sample=samples.index),
        expand('plots/png/{sample}_knee_plot.png', sample=samples.index),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv',
        'reports/combined.html',
        'plots/umi_per_gene_distribution.pdf',
        'plots/png/umi_per_gene_distribution.png'

rule meta:
    input:
        expand('{annotation_prefix}.refFlat',annotation_prefix=annotation_prefix),
        expand('{annotation_prefix}.reduced.gtf',annotation_prefix=annotation_prefix),
        expand('{reference_prefix}.dict',reference_prefix=reference_prefix),
        expand('{reference_prefix}.rRNA.intervals',reference_prefix=reference_prefix),
        expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths)

rule qc:
    input:
        expand('logs/{sample}_R1_fastqc.html', sample=samples.index),
        'reports/fastqc.html',
        'reports/fastqc_data/multiqc_general_stats.txt'

rule filter:
    input:
        expand('data/{sample}_filtered.fastq.gz', sample=samples.index),
        expand('plots/{sample}_polya_trimmed.pdf', sample=samples.index),
        expand('plots/{sample}_start_trim.pdf', sample=samples.index),
        expand('plots/{sample}_CELL_dropped.pdf', sample=samples.index),
        expand('plots/{sample}_UMI_dropped.pdf', sample=samples.index),
        'reports/filter.html',
        'plots/BC_drop.pdf',
        'plots/png/BC_drop.png'

rule map:
    input:  
        expand('data/{sample}_final.bam', sample=samples.index),
        expand('logs/{sample}_hist_out_cell.txt', sample=samples.index),
        expand('plots/{sample}_knee_plot.pdf', sample=samples.index),
        expand('plots/png/{sample}_knee_plot.png', sample=samples.index),
        'reports/star.html'

rule extract:
    input:
        expand('logs/{sample}_umi_per_gene.tsv', sample=samples.index),
        'summary/umi_expression_matrix.tsv',
        'summary/counts_expression_matrix.tsv',
        expand('plots/{sample}_rna_metrics.pdf', sample=samples.index),
        expand('plots/png/{sample}_rna_metrics.png', sample=samples.index),
        'plots/umi_per_gene_distribution.pdf',
        'plots/png/umi_per_gene_distribution.png'
        

rule split_species:
    input:
        expand('summary/{sample}_{species}_barcodes.csv', sample=samples.index, species=config['META']['species']),
        expand('plots/{sample}_species_plot_genes.pdf', sample=samples.index),
        expand('plots/{sample}_species_plot_transcripts.pdf', sample=samples.index),
        expand('plots/png/{sample}_species_plot_genes.png', sample=samples.index),
        expand('plots/png/{sample}_species_plot_transcripts.png', sample=samples.index)
        
        

rule extract_species:
    input:
        expand('summary/{sample}_{species}_umi_expression_matrix.txt', sample=samples.index, species=config['META']['species']),
        expand('summary/{sample}_{species}_counts_expression_matrix.txt', sample=samples.index, species=config['META']['species']),
        expand('logs/{sample}_{species}_umi_per_gene.tsv', sample=samples.index, species=config['META']['species']),
        expand('summary/Experiment_{species}_counts_expression_matrix.tsv', species=config['META']['species']),
        expand('summary/Experiment_{species}_umi_expression_matrix.tsv', species=config['META']['species']),
        expand('plots/species_{sample}_{species}_rna_metrics.pdf', sample=samples.index, species=config['META']['species']),
        expand('plots/png/species_{sample}_{species}_rna_metrics.png', sample=samples.index, species=config['META']['species'])


rule multiqc_all:
    input: 
        fastqc = expand('logs/{sample}_R1_fastqc.html',sample=samples.index),
        trimmomatic=expand('logs/{sample}_trimlog.txt',sample=samples.index),
        star=expand('logs/{sample}.Log.final.out',sample=samples.index)
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