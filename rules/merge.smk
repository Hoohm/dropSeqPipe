
localrules:
	merge_long,
	violine_plots

rule merge_long:
    input:
        expand('{results_dir}/samples/{sample}/{{type}}/expression.long', sample=samples.index, results_dir=results_dir)
    output:
        mtx='{results_dir}/summary/{type}/expression.mtx',
        barcodes='{results_dir}/summary/{type}/barcodes.tsv',
        features='{results_dir}/summary/{type}/features.tsv',
    params:
        samples=lambda wildcards: samples.index
    script:
        "../scripts/convert_mtx.py"


rule violine_plots:
    input:
        UMIs='{results_dir}/summary/umi/expression.mtx',
        counts='{results_dir}/summary/read/expression.mtx',
        design='samples.csv'
    conda: '../envs/plots_ext.yaml'
    output:
        pdf_violine='{results_dir}/plots/violinplots_comparison_UMI.pdf',
        pdf_umivscounts='{results_dir}/plots/UMI_vs_counts.pdf',
        pdf_umi_vs_gene='{results_dir}/plots/UMI_vs_gene.pdf',
        pdf_count_vs_gene='{results_dir}/plots/Count_vs_gene.pdf',
        R_objects='{results_dir}/summary/R_Seurat_objects.rdata'
    script:
        '../scripts/plot_violine.R'
