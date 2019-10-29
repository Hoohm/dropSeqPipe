

ruleorder: extend_barcode_whitelist > extend_barcode_top
ruleorder: extend_barcode_whitelist > get_cell_whitelist


localrules:
    get_cell_whitelist,
    extend_barcode_top

rule extend_barcode_whitelist:
    input:
        whitelist=barcode_whitelist
    output:
        barcodes='{results_dir}/samples/{sample}/barcodes.csv',
        barcode_ref='{results_dir}/samples/{sample}/barcode_ref.pkl',
        barcode_ext_ref='{results_dir}/samples/{sample}/barcode_ext_ref.pkl',
        barcode_mapping='{results_dir}/samples/{sample}/empty_barcode_mapping.pkl'
    script:
        '../scripts/generate_extended_ref.py'

rule get_top_barcodes:
    input:
        '{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz'
    output:
        '{results_dir}/samples/{sample}/top_barcodes.csv'
    conda: '../envs/umi_tools.yaml'
    params:
        cell_barcode_length=(config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1),
        umi_barcode_length=(config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'] + 1),
        num_cells=lambda wildcards: round(int(samples.loc[wildcards.sample,'expected_cells'])*1.2),
    shell:
        """umi_tools whitelist\
        --stdin {input}\
        --bc-pattern='(?P<cell_1>.{{{params.cell_barcode_length}}})(?P<umi_1>.{{{params.umi_barcode_length}}})'\
        --extract-method=regex\
        --set-cell-number={params.num_cells}\
        --log2stderr > {output}"""

rule get_cell_whitelist:
    input:
        '{results_dir}/samples/{sample}/top_barcodes.csv'
    output:
        '{results_dir}/samples/{sample}/barcodes.csv'
    shell:
        """cat {input} | cut -f 1 > {output}"""


rule extend_barcode_top:
    input:
        whitelist='{results_dir}/samples/{sample}/top_barcodes.csv'
    output:
        barcode_ref='{results_dir}/samples/{sample}/barcode_ref.pkl',
        barcode_ext_ref='{results_dir}/samples/{sample}/barcode_ext_ref.pkl',
        barcode_mapping='{results_dir}/samples/{sample}/empty_barcode_mapping.pkl'
    script:
        '../scripts/umi_tools_extended_ref.py'


rule repair_barcodes:
    input:
        bam='{results_dir}/samples/{sample}/Aligned.merged.bam',
        barcode_ref='{results_dir}/samples/{sample}/barcode_ref.pkl',
        barcode_ext_ref='{results_dir}/samples/{sample}/barcode_ext_ref.pkl',
        barcode_mapping='{results_dir}/samples/{sample}/empty_barcode_mapping.pkl'
    conda: '../envs/merge_bam.yaml'
    output:
        bam=temp('{results_dir}/samples/{sample}/Aligned.repaired.bam'),
        barcode_mapping_counts='{results_dir}/samples/{sample}/barcode_mapping_counts.pkl'
    script:
        '../scripts/repair_barcodes.py'