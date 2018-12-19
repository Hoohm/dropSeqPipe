"""Extract species specific expression to prepare the species plot."""


#Which rules will be run on the host computer and not sent to nodes
localrules: plot_barnyard

rule split_bam_species:
    input:
        '{results_dir}/samples/{sample}/final.bam'
    output:
        '{results_dir}/samples/{sample}/{species}/unfiltered.bam'
    params:
        species=lambda wildcards: wildcards.species,
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && FilterBam -m {params.memory}\
        REF_SOFT_MATCHED_RETAINED={params.species}\
        INPUT={input}\
        OUTPUT={output}"""


rule extract_all_umi_expression:
    input: 
        data='{results_dir}/samples/{sample}/{species}/unfiltered.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/barcodes.csv'
    output:
        umi_matrix=temp('{results_dir}/samples/{sample}/{species}/unfiltered_umi_expression_matrix.tsv'),
        summary='{results_dir}/samples/{sample}/{species}/dge.summary.txt'
    params:
        count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
        cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory'],
        locus_list=','.join(config['EXTRACTION']['LOCUS'])
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
        I={input.data}\
        O={output.umi_matrix}\
        SUMMARY={output.summary}\
        EDIT_DISTANCE={params.cellBarcodeEditDistance}\
        CELL_BC_FILE={input.barcode_whitelist}\
        LOCUS_FUNCTION_LIST={{{params.locus_list}}}\
        MIN_BC_READ_THRESHOLD={params.count_per_umi}"""


rule plot_barnyard:
    input:
        expand('{results_dir}/samples/{{sample}}/{species}/dge.summary.txt',species=config['META']['species'], results_dir=results_dir)
    output: 
        genes_pdf='{results_dir}/plots/barnyard/{sample}_genes.pdf',
        transcripts_pdf='{results_dir}/plots/barnyard/{sample}_transcripts.pdf',
        barcodes_species=expand('{{results_dir}}/samples/{{sample}}/{species}/barcodes.csv', species=species_list)
    params:
        expected_cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells'])
    script: 
        '../scripts/plot_species_plot.R'