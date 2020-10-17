localrules:
    violine_plots,
    summary_stats

# rule merge_long:
#     input:
#         expand('{results_dir}/samples/{sample}/{{type}}/expression.long', sample=samples.index, results_dir=results_dir)
#     output:
#         mtx='{results_dir}/summary/{type}/matrix.mtx',
#         barcodes='{results_dir}/summary/{type}/barcodes.tsv',
#         features='{results_dir}/summary/{type}/features.tsv',
#     params:
#         samples=lambda wildcards: samples.index
#     conda: '../envs/merge_long.yaml'
#     script:
#         "../scripts/convert_mtx.py"

rule violine_plots:
    input:
        umi_mtx=expand('{results_dir}/summary/{sample}/umi/matrix.mtx.gz', results_dir=results_dir, sample=samples.index),
        design='samples.csv'
    conda: '../envs/r.yaml'
    output:
        pdf_violine='{results_dir}/plots/violinplots_comparison_UMI.pdf',
        pdf_umivscounts='{results_dir}/plots/UMI_vs_counts.pdf',
        pdf_umi_vs_gene='{results_dir}/plots/UMI_vs_gene.pdf',
        pdf_count_vs_gene='{results_dir}/plots/Count_vs_gene.pdf',
        R_objects='{results_dir}/summary/R_Seurat_objects.rdata'
    script:
        '../scripts/plot_violine.R'

rule summary_stats:
    input:
        R_objects='{results_dir}/summary/R_Seurat_objects.rdata',
        R2qc=expand('{results_dir}/logs/cutadapt/{sample}_R2.qc.txt', sample=samples.index, results_dir=results_dir),
        hist_cell=expand('{results_dir}/logs/dropseq_tools/{sample}_hist_out_cell.txt', sample=samples.index, results_dir=results_dir)
    conda: '../envs/r.yaml'
    output:
        stats_pre='{results_dir}/summary/barcode_stats_pre_filter.csv',
        stats_post='{results_dir}/summary/barcode_stats_post_filter.csv',
    params:
        sample_names=lambda wildcards: samples.index,
        batches=lambda wildcards: samples.loc[samples.index, 'batch']
    script:
        '../scripts/create_summary_stats.R'
