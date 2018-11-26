import pandas as pd
import os
import re
import glob


# Load configuration file
configfile: "config.yaml"


# Get sample names from samples.csv
samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)
types=['reads','umi']
# Get read_lengths from samples.csv
read_lengths = list(samples.loc[:,'read_length'])

# Constraint sample names wildcards
wildcard_constraints:
    sample="({})".format("|".join(samples.index))


results_dir = config['LOCAL']['results']
raw_data_dir = config['LOCAL']['raw_data']



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
        expand(
        ['{annotation_prefix}.refFlat',
        '{annotation_prefix}.reduced.gtf',
        '{reference_prefix}.dict',
        '{reference_prefix}.rRNA.intervals',
        '{star_index_prefix}_{read_length}/SA',
        #qc
        '{results_dir}reports/fastqc_reads.html',
        '{results_dir}reports/fastqc_barcodes.html',
        #filter
        '{results_dir}plots/adapter_content.pdf',
        '{results_dir}reports/barcode_filtering.html',
        '{results_dir}reports/RNA_filtering.html',
        #mapping
        '{results_dir}plots/knee_plots/{sample}_knee_plot.pdf',
        '{results_dir}reports/star.html',
        '{results_dir}plots/yield.pdf',
        #extract
        '{results_dir}plots/rna_metrics/{sample}_rna_metrics.pdf',
        '{results_dir}summary/{type}/expression.mtx',
        '{results_dir}samples/{sample}/{type}/expression.mtx',
        #merge
        '{results_dir}plots/UMI_vs_counts.pdf',
        '{results_dir}plots/UMI_vs_gene.pdf',
        '{results_dir}plots/Count_vs_gene.pdf',
        '{results_dir}summary/R_Seurat_objects.rdata',
        '{results_dir}plots/violinplots_comparison_UMI.pdf'],
            annotation_prefix=annotation_prefix,
            reference_prefix=reference_prefix,
            star_index_prefix=star_index_prefix,
            read_length=read_lengths,
            sample=samples.index,
            type=types,
            results_dir=results_dir)
        
rule test:
    input:
        expand("samples/{sample}/trimmmed_R1.fastq.gz", sample=samples.index)
        
rule meta:
    input:
        expand(
        ['{annotation_prefix}.refFlat',
        '{annotation_prefix}.reduced.gtf',
        '{reference_prefix}.dict',
        '{reference_prefix}.rRNA.intervals',
        '{star_index_prefix}_{read_length}/SA'],
        star_index_prefix=star_index_prefix,
        read_length=read_lengths,
        reference_prefix=reference_prefix,
        annotation_prefix=annotation_prefix)

rule qc:
    input:
        expand(
        ['{results_dir}reports/fastqc_reads.html',
        '{results_dir}reports/fastqc_barcodes.html'],
        results_dir=results_dir)

rule filter:
    input:
        expand(
        ['{results_dir}plots/adapter_content.pdf',
        '{results_dir}reports/barcode_filtering.html',
        '{results_dir}reports/RNA_filtering.html',
        '{results_dir}samples/{sample}/trimmmed_repaired_R1.fastq.gz'],
        results_dir=results_dir,
        sample=samples.index)
        
rule map:
    input:  
        expand(
        ['{results_dir}plots/knee_plots/{sample}_knee_plot.pdf',
        '{results_dir}reports/star.html',
        '{results_dir}plots/yield.pdf',
        '{results_dir}samples/{sample}/final.bam'],
        sample=samples.index,
        results_dir=results_dir)

rule extract:
    input:
        expand(
        ['{results_dir}plots/rna_metrics/{sample}_rna_metrics.pdf',
        '{results_dir}summary/{type}/expression.mtx',
        '{results_dir}samples/{sample}/{type}/expression.mtx'],
        results_dir=results_dir,
        sample=samples.index,
        type=types)

rule split_species:
    input:
        expand(
        ['{results_dir}samples/{sample}/{species}/barcodes.csv',
        '{results_dir}plots/barnyard/{sample}_genes.pdf',
        '{results_dir}plots/barnyard/{sample}_transcripts.pdf',
        '{results_dir}samples/{sample}/{species}/unfiltered.bam'],
        sample=samples.index,
        species=config['META']['species'],
        results_dir=results_dir)


rule extract_species:
    input:
        expand(
        ['{results_dir}samples/{sample}/{species}/umi_expression_matrix.txt',
        '{results_dir}samples/{sample}/{species}/counts_expression_matrix.txt',
        '{results_dir}summary/Experiment_{species}_counts_expression_matrix.tsv',
        '{results_dir}summary/Experiment_{species}_umi_expression_matrix.tsv',
        '{results_dir}plots/rna_metrics/{sample}_{species}_rna_metrics.pdf'],
        sample=samples.index,
        species=config['META']['species'],
        results_dir=results_dir)
        
rule merge:
    input:
        #merge
        expand(
        ['{results_dir}plots/UMI_vs_counts.pdf',
        '{results_dir}plots/UMI_vs_gene.pdf',
        '{results_dir}plots/Count_vs_gene.pdf',
        '{results_dir}summary/R_Seurat_objects.rdata',
        '{results_dir}plots/violinplots_comparison_UMI.pdf'],
        results_dir=results_dir)
        


include: "rules/generate_meta.smk"
include: "rules/fastqc.smk"
include: "rules/filter.smk"
include: "rules/cell_barcodes.smk"
include: "rules/map.smk"
include: "rules/extract_expression_single.smk"
include: "rules/split_species.smk"
include: "rules/extract_expression_species.smk"
include: "rules/merge.smk"
