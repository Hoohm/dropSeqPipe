localrules:
    violine_plots,
    summary_stats,
    create_seurat_object


rule create_seurat_and_sce_object:
    input:
        umi_mtx=expand('{results_dir}/samples/{sample}/umi/matrix.mtx.gz', results_dir=results_dir, sample=samples.index),
        rna_metrics=expand('{results_dir}/logs/dropseq_tools/{sample}_rna_metrics.txt', results_dir=results_dir, sample=samples.index),
    output:
        seurat_object='{results_dir}/summary/Seurat_object.rds',
        SCE_object='{results_dir}/summary/SCE_object.rds'
    params:
        design='samples.csv',
        samples = samples.index
    conda:
        '../envs/r.yaml'
    log:
        stdout='{results_dir}/logs/merge/create_seurat_object.Rout',
    script:
        '../scripts/create_seurat_object.R'

rule violine_plots:
    input:
        seurat_object='{results_dir}/summary/Seurat_object.rds',
        design='samples.csv'
    conda: '../envs/r.yaml'
    output:
        pdf_violine='{results_dir}/plots/violinplots_comparison_UMI.pdf',
        pdf_umivscounts='{results_dir}/plots/UMI_vs_counts.pdf',
        pdf_umi_vs_gene='{results_dir}/plots/UMI_vs_gene.pdf',
        pdf_count_vs_gene='{results_dir}/plots/Count_vs_gene.pdf',
    script:
        '../scripts/plot_violine.R'

# rule summary_stats:
#     input:
#         R_objects='{results_dir}/summary/Seurat_object.rds',
#         R2qc=expand('{results_dir}/logs/cutadapt/{sample}_R2.qc.txt', sample=samples.index, results_dir=results_dir),
#         hist_cell=expand('{results_dir}/logs/dropseq_tools/{sample}_hist_out_cell.txt', sample=samples.index, results_dir=results_dir)
#     conda: '../envs/r.yaml'
#     output:
#         stats_pre='{results_dir}/summary/barcode_stats_pre_filter.csv',
#         stats_post='{results_dir}/summary/barcode_stats_post_filter.csv',
#     params:
#         sample_names=lambda wildcards: samples.index,
#         batches=lambda wildcards: samples.loc[samples.index, 'batch']
#     script:
#         '../scripts/create_summary_stats.R'
