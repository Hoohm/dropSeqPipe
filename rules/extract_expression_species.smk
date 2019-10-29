"""Extract expression fof mixed species"""

#Which rules will be run on the host computer and not sent to nodes
localrules:
    plot_rna_metrics_species,
    convert_long_to_mtx_species

rule extract_umi_expression_species:
    input:
        data='{results_dir}/samples/{sample}/{species}/unfiltered.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/{species}/barcodes.csv'
    output:
        dense=temp('{results_dir}/samples/{sample}/{species}/umi/expression.txt'),
        long=temp('{results_dir}/samples/{sample}/{species}/umi/expression.long')
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


rule extract_reads_expression_species:
    input:
        data='{results_dir}/samples/{sample}/{species}/unfiltered.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/{species}/barcodes.csv'
    params:
        count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
        num_cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),
        umiBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory'],
        locus_list=','.join(config['EXTRACTION']['LOCUS']),
        strand_strategy=config['EXTRACTION']['strand-strategy']
    output:
        dense=temp('{results_dir}/samples/{sample}/{species}/read/expression.txt'),
        long=temp('{results_dir}/samples/{sample}/{species}/read/expression.long')
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

rule convert_long_to_mtx_species:
    input:
        '{results_dir}/samples/{sample}/{species}/{type}/expression.long'
    output:
        barcodes='{results_dir}/samples/{sample}/{species}/{type}/barcodes.tsv',
        features='{results_dir}/samples/{sample}/{species}/{type}/genes.tsv',
        mtx='{results_dir}/samples/{sample}/{species}/{type}/matrix.mtx'
    params:
        samples=lambda wildcards: wildcards.sample
    script:
        "../scripts/convert_mtx.py"

rule compress_mtx_species:
    input: 
        barcodes='{results_dir}/samples/{sample}/{species}/{type}/barcodes.tsv',
        features='{results_dir}/samples/{sample}/{species}/{type}/genes.tsv',
        mtx='{results_dir}/samples/{sample}/{species}/{type}/matrix.mtx'
    output:
        barcodes='{results_dir}/samples/{sample}/{species}/{type}/barcodes.tsv.gz',
        features='{results_dir}/samples/{sample}/{species}/{type}/genes.tsv.gz',
        mtx='{results_dir}/samples/{sample}/{species}/{type}/matrix.mtx.gz'
    conda: '../envs/pigz.yaml'
    threads: 3
    shell:
        """pigz -p {threads} {input.barcodes} {input.features} {input.mtx}"""

rule SingleCellRnaSeqMetricsCollector_species:
    input:
        data='{results_dir}/samples/{sample}/{species}/unfiltered.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/{species}/barcodes.csv',
        refFlat=expand("{ref_path}/{species}_{build}_{release}/curated_annotation.refFlat",
            ref_path=ref_path,
            release=release,
            species=species,
            build=build),
        rRNA_intervals=expand("{ref_path}/{species}_{build}_{release}/annotation.rRNA.intervals",
            ref_path=ref_path,
            release=release,
            build=build,
            species=species)
    params:
        cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),        
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        '{results_dir}/logs/dropseq_tools/{sample}/{species}/rna_metrics.txt'
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && SingleCellRnaSeqMetricsCollector -m {params.memory}\
        INPUT={input.data}\
        OUTPUT={output}\
        ANNOTATIONS_FILE={input.refFlat}\
        CELL_BC_FILE={input.barcode_whitelist}\
        RIBOSOMAL_INTERVALS={input.rRNA_intervals}
        """
rule plot_rna_metrics_species:
    input:
        rna_metrics='{results_dir}/logs/dropseq_tools/{sample}/{species}/rna_metrics.txt',
        barcode='{results_dir}/samples/{sample}/{species}/barcodes.csv'
    conda: '../envs/r.yaml'
    output:
        pdf='{results_dir}/plots/rna_metrics/{sample}_{species}_rna_metrics.pdf'
    script:
        '../scripts/plot_rna_metrics.R'
