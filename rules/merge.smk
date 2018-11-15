

rule convert_long_to_mtx:
    input:
        'summary/{sample}_{type}_expression.long'
    output:
        barcodes='summary/{sample}/{type}/barcodes.tsv',
        features='summary/{sample}/{type}/features.tsv',
        mtx='summary/{sample}/{type}/expression.mtx'
    params:
        samples=lambda wildcards: wildcards.sample
    script:
        "../scripts/convert_mtx.py"

rule merge_long:
    input:
        expand('summary/{sample}_{{type}}_expression.long', sample=samples.index)
    output:
        mtx='summary/experiment/{type}/expression.mtx',
        barcodes='summary/experiment/{type}/barcodes.tsv',
        features='summary/experiment/{type}/features.tsv',
    params:
        samples=lambda wildcards: samples.index
    script:
        "../scripts/convert_mtx.py"


rule violine_plots:
    input:
        UMIs='summary/experiment/umi/expression.mtx',
        counts='summary/experiment/reads/expression.mtx',
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