"""Extract expression fof single species"""

#Which rules will be run on the host computer and not sent to nodes
localrules: plot_rna_metrics

rule extract_umi_expression:
    input:
        data='data/{sample}/final.bam',
        barcode_whitelist='data/{sample}/barcodes.csv'
    output:
        sparse='data/{sample}/umi/expression.long',
        dense=temp('data/{sample}/umi/expression.tsv')
    params:
        count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
        num_cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),
        umiBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
        I={input.data}\
        O={output.dense}\
        EDIT_DISTANCE={params.umiBarcodeEditDistance}\
        OUTPUT_LONG_FORMAT={output.sparse}\
        STRAND_STRATEGY=SENSE\
        OUTPUT_READS_INSTEAD=false\
        MIN_BC_READ_THRESHOLD={params.count_per_umi}\
        CELL_BC_FILE={input.barcode_whitelist}"""

rule extract_reads_expression:
    input:
        data='data/{sample}/final.bam',
        barcode_whitelist='data/{sample}/barcodes.csv'
    output:
        sparse='data/{sample}/reads/expression.long',
        dense=temp('data/{sample}/reads/expression.tsv')
    params:
        count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
        num_cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),
        umiBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
        I={input.data}\
        O={output.dense}\
        EDIT_DISTANCE={params.umiBarcodeEditDistance}\
        OUTPUT_LONG_FORMAT={output.sparse}\
        STRAND_STRATEGY=SENSE\
        OUTPUT_READS_INSTEAD=true\
        MIN_BC_READ_THRESHOLD={params.count_per_umi}\
        CELL_BC_FILE={input.barcode_whitelist}"""


rule SingleCellRnaSeqMetricsCollector:
    input:
        data='data/{sample}/final.bam',
        barcode_whitelist='data/{sample}/barcodes.csv',
        refFlat="{}.refFlat".format(annotation_prefix),
        rRNA_intervals="{}.rRNA.intervals".format(reference_prefix)
    params:     
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory']
    output:
        rna_metrics='logs/dropseq_tools/{sample}_rna_metrics.txt',
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
        rna_metrics='logs/dropseq_tools/{sample}_rna_metrics.txt',
        barcodes='data/{sample}/barcodes.csv'
    conda: '../envs/plots.yaml'
    output:
        pdf='plots/rna_metrics/{sample}_rna_metrics.pdf'
    script:
        '../scripts/plot_rna_metrics.R'