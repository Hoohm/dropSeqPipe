
localrules: convert_long_to_mtx, merge_long, violine_plots

rule convert_long_to_mtx:
    input:
        'data/{sample}/{type}/expression.long'
    output:
        barcodes='data/{sample}/{type}/barcodes.tsv',
        features='data/{sample}/{type}/features.tsv',
        mtx='data/{sample}/{type}/expression.mtx'
    params:
        samples=lambda wildcards: wildcards.sample
    script:
        "../scripts/convert_mtx.py"

rule merge_long:
    input:
        expand('data/{sample}/{{type}}/expression.long', sample=samples.index)
    output:
        mtx='summary/{type}/expression.mtx',
        barcodes='summary/{type}/barcodes.tsv',
        features='summary/{type}/features.tsv',
    params:
        samples=lambda wildcards: samples.index
    script:
        "../scripts/convert_mtx.py"


rule violine_plots:
    input:
        UMIs='summary/umi/expression.mtx',
        counts='summary/reads/expression.mtx',
        design='samples.csv'
    conda: '../envs/plots_ext.yaml'
    output:
        pdf_violine='plots/violinplots_comparison_UMI.pdf',
        pdf_umivscounts='plots/UMI_vs_counts.pdf',
        pdf_umi_vs_gene='plots/UMI_vs_gene.pdf',
        pdf_count_vs_gene='plots/Count_vs_gene.pdf',
        R_objects='summary/R_Seurat_objects.rdata'
    script:
        '../scripts/plot_violine.R'