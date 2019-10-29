"""Extract expression fof single species"""

#Which rules will be run on the host computer and not sent to nodes
localrules:
    plot_rna_metrics,
    convert_long_to_mtx,
    compress_mtx

rule extract_umi_expression:
    input:
        data='{results_dir}/samples/{sample}/final.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/barcodes.csv'
    output:
        long='{results_dir}/samples/{sample}/umi/expression.long',
        dense=temp('{results_dir}/samples/{sample}/umi/expression.tsv')
    params:
        count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
        num_cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),
        umiBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory'],
        locus_list=','.join(config['EXTRACTION']['LOCUS']),
        strand_strategy=config['EXTRACTION']['strand-strategy']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
        I={input.data}\
        O={output.dense}\
        EDIT_DISTANCE={params.umiBarcodeEditDistance}\
        OUTPUT_LONG_FORMAT={output.long}\
        STRAND_STRATEGY={params.strand_strategy}\
        OUTPUT_READS_INSTEAD=false\
        LOCUS_FUNCTION_LIST={{{params.locus_list}}}\
        MIN_BC_READ_THRESHOLD={params.count_per_umi}\
        CELL_BC_FILE={input.barcode_whitelist}"""

rule extract_reads_expression:
    input:
        data='{results_dir}/samples/{sample}/final.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/barcodes.csv'
    output:
        long=temp('{results_dir}/samples/{sample}/read/expression.long'),
        dense=temp('{results_dir}/samples/{sample}/read/expression.tsv')
    params:
        count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
        num_cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),
        umiBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory'],
        locus_list=','.join(config['EXTRACTION']['LOCUS']),
        strand_strategy=config['EXTRACTION']['strand-strategy']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
        I={input.data}\
        O={output.dense}\
        EDIT_DISTANCE={params.umiBarcodeEditDistance}\
        OUTPUT_LONG_FORMAT={output.long}\
        STRAND_STRATEGY={params.strand_strategy}\
        OUTPUT_READS_INSTEAD=true\
        LOCUS_FUNCTION_LIST={{{params.locus_list}}}\
        MIN_BC_READ_THRESHOLD={params.count_per_umi}\
        CELL_BC_FILE={input.barcode_whitelist}"""


rule SingleCellRnaSeqMetricsCollector:
    input:
        data='{results_dir}/samples/{sample}/final.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/barcodes.csv',
        refFlat=expand("{ref_path}/{species}_{build}_{release}/curated_annotation.refFlat",
            ref_path=config['META']['reference-directory'],
            species=species,
            release=release,
            build=build),
        rRNA_intervals=expand("{ref_path}/{species}_{build}_{release}/annotation.rRNA.intervals",
            ref_path=config['META']['reference-directory'],
            species=species,
            release=release,
            build=build)
    params:     
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory']
    output:
        rna_metrics='{results_dir}/logs/dropseq_tools/{sample}_rna_metrics.txt',
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && SingleCellRnaSeqMetricsCollector -m {params.memory}\
        INPUT={input.data}\
        OUTPUT={output}\
        ANNOTATIONS_FILE={input.refFlat}\
        CELL_BC_FILE={input.barcode_whitelist}\
        RIBOSOMAL_INTERVALS={input.rRNA_intervals}
        """

rule plot_rna_metrics:
    input:
        rna_metrics='{results_dir}/logs/dropseq_tools/{sample}_rna_metrics.txt',
        barcodes='{results_dir}/samples/{sample}/barcodes.csv'
    conda: '../envs/r.yaml'
    output:
        pdf='{results_dir}/plots/rna_metrics/{sample}_rna_metrics.pdf'
    script:
        '../scripts/plot_rna_metrics.R'


rule convert_long_to_mtx:
    input:
        '{results_dir}/samples/{sample}/{type}/expression.long'
    output:
        barcodes='{results_dir}/samples/{sample}/{type}/barcodes.tsv',
        features='{results_dir}/samples/{sample}/{type}/genes.tsv',
        mtx='{results_dir}/samples/{sample}/{type}/matrix.mtx'
    params:
        samples=lambda wildcards: wildcards.sample
    script:
        "../scripts/convert_mtx.py"

rule compress_mtx:
    input: 
        barcodes='{results_dir}/samples/{sample}/{type}/barcodes.tsv',
        features='{results_dir}/samples/{sample}/{type}/genes.tsv',
        mtx='{results_dir}/samples/{sample}/{type}/matrix.mtx'
    output:
        barcodes='{results_dir}/samples/{sample}/{type}/barcodes.tsv.gz',
        features='{results_dir}/samples/{sample}/{type}/genes.tsv.gz',
        mtx='{results_dir}/samples/{sample}/{type}/matrix.mtx.gz'
    conda: '../envs/pigz.yaml'
    threads: 3
    shell:
        """pigz -p {threads} {input.barcodes} {input.features} {input.mtx}"""