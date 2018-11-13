
rule violine_plots:
    input:
        UMIs='summary/umi_expression_matrix.tsv',
        counts='summary/counts_expression_matrix.tsv',
        design='samples.csv'
#   params:
#       cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
#       edit_distance=config['EXTRACTION']['UMI-edit-distance']
    conda: '../envs/plots_ext.yaml'
    output:
        pdf_violine='plots/violinplots_comparison_UMI.pdf',
#       html_umivscounts='plots/UMI_vs_counts.html',
        pdf_umivscounts='plots/UMI_vs_counts.pdf',
#       html_umi_vs_gene='plots/UMI_vs_gene.html',
        pdf_umi_vs_gene='plots/UMI_vs_gene.pdf',
#       html_count_vs_gene='plots/Count_vs_gene.html',
        pdf_count_vs_gene='plots/Count_vs_gene.pdf',
        R_objects='summary/R_Seurat_objects.rdata'
    script:
        '../scripts/plot_violine.R'

rule convert_long_to_mtx:
    input:
        long_file='summary/{sample}_{type}_expression.long'
    output:
        barcodes='summary/{sample}/{type}/barcodes.csv',
        features='summary/{sample}/{type}/features.csv',
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
        barcodes='summary/experiment/{type}/barcodes.csv',
        features='summary/experiment/{type}/features.csv'
    params:
        samples=lambda wildcards: samples.index
    script:
        "../scripts/convert_mtx.py"

rule merge_umi:
    input:
        expand('summary/{sample}_umi_expression_matrix.tsv', sample=samples.index)
    params:
        sample_names=lambda wildcards: samples.index
    conda: '../envs/merge.yaml'
    output:
        'summary/umi_expression_matrix.tsv'
    script:
        "../scripts/merge_counts_single.R"

rule merge_counts:
    input:
        expand('summary/{sample}_counts_expression_matrix.tsv', sample=samples.index)
    params:
        sample_names=lambda wildcards: samples.index
    conda: '../envs/merge.yaml'
    output:
        'summary/counts_expression_matrix.tsv'
    script:
        "../scripts/merge_counts_single.R"
