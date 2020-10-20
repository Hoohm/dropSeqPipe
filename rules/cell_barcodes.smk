


ruleorder: copy_whitelist > get_cell_whitelist


localrules:
    get_cell_whitelist,
    copy_whitelist,
    filter_top_barcodes


rule get_top_barcodes:
    input:
        '{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz'
    output:
        '{results_dir}/samples/{sample}/umitools_barcode_whitelist.csv'
    log:
        '{results_dir}/logs/{sample}/umitools.log'
    conda: '../envs/umi_tools.yaml'
    params:
        cell_barcode_length=(config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1),
        umi_barcode_length=(config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'] + 1),
        num_cells=lambda wildcards: round(int(samples.loc[wildcards.sample,'expected_cells'])),
    shell:
        """umi_tools whitelist\
        --stdin {input}\
        --bc-pattern='(?P<cell_1>.{{{params.cell_barcode_length}}})(?P<umi_1>.{{{params.umi_barcode_length}}})'\
        --extract-method=regex\
        --set-cell-number={params.num_cells}\
        --plot-prefix={wildcards.results_dir}/plots/umitools_{wildcards.sample}\
        --log2stderr\
        --log={log}\
        --stdout={output}
        """

rule get_cell_whitelist:
    input:
        '{results_dir}/samples/{sample}/umitools_barcode_whitelist.csv'
    output:
        '{results_dir}/samples/{sample}/barcodes.csv'
    shell:
        """cat {input} | cut -f 1 > {output}"""

rule copy_whitelist:
    input:
        whitelist = barcode_whitelist
    output:
        '{results_dir}/samples/{sample}/barcodes.csv'
    shell:
        """cp {input.whitelist} {output}"""


# rule bam_hist:
#     input:
#         '{results_dir}/samples/{sample}/final.bam'
#     params:
#         memory=config['LOCAL']['memory'],
#         temp_directory=config['LOCAL']['temp-directory']
#     output:
#         '{results_dir}/logs/dropseq_tools/{sample}_hist_out_cell.txt'
#     conda: '../envs/dropseq_tools.yaml'
#     shell:
#         """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && BamTagHistogram -m {params.memory}\
#         TAG=CB\
#         I={input}\
#         READ_MQ=10\
#         O={output}
#         """

#rule filter_top_barcodes:
#    input:
#        whitelist='{results_dir}/samples/{sample}/barcodes.csv',
#        top_barcodes='{results_dir}/logs/dropseq_tools/{sample}_hist_out_cell.txt'
#    output:
#        filtered_barcodes='{results_dir}/samples/{sample}/filtered_barcodes.csv'
#    log:
#        filter_log='{results_dir}/logs/custom/{sample}_unrecognized_barcodes.csv'
#    params:
#        num_cells=lambda wildcards: round(int(samples.loc[wildcards.sample,'expected_cells']))
#    # enforce conde env to ensure consistent python version for python scripts
#    conda: '../envs/cutadapt.yaml'
#    script:
#        '../scripts/whitelist_check.py'
