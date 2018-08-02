"""Extract species specific expression to prepare the species plot."""

ruleorder: extract_all_umi_expression_whitelist_species > extract_all_umi_expression_species

#Which rules will be run on the host computer and not sent to nodes
localrules: plot_barnyard

rule split_bam_species:
	input:
		'data/{sample}_final.bam'
	output:
		'data/{species}/{sample}_unfiltered.bam'
	params:
		species=lambda wildcards: wildcards.species,
		dropseq_wrapper='../scripts/drop-seq-tools-wrapper.sh',
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""{params.dropseq_wrapper} -t {params.temp_directory} -m {params.memory} -p FilterBAM\
		REF_SOFT_MATCHED_RETAINED={params.species}\
		INPUT={input}\
		OUTPUT={output}"""


rule extract_all_umi_expression_species:
	input:
		'data/{species}/{sample}_unfiltered.bam'
	output:
		umi_matrix=temp('summary/{species}/{sample}_unfiltered_umi_expression_matrix.tsv'),
		summary='summary/{species}/{sample}_dge.summary.txt'
	params:
		count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
		num_cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
		dropseq_wrapper='../scripts/drop-seq-tools-wrapper.sh',
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""{params.dropseq_wrapper} -t {params.temp_directory} -m {params.memory} -p DigitalExpression\
		I={input}\
		O={output.umi_matrix}\
		SUMMARY={output.summary}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		NUM_CORE_BARCODES={params.num_cells}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""

rule extract_all_umi_expression_whitelist_species:
	input: 
		data='data/{species}/{sample}_unfiltered.bam',
		barcode_whitelist='barcodes.csv'
	output:
		umi_matrix=temp('summary/{species}/{sample}_unfiltered_umi_expression_matrix.tsv'),
		summary='summary/{species}/{sample}_dge.summary.txt'
	params:
		count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
		dropseq_wrapper='../scripts/drop-seq-tools-wrapper.sh',
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""{params.dropseq_wrapper} -t {params.temp_directory} -m {params.memory} -p DigitalExpression\
		I={input.data}\
		O={output.umi_matrix}\
		SUMMARY={output.summary}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		CELL_BC_FILE={input.barcode_whitelist}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""


rule plot_barnyard:
	input:
		expand('summary/{species}/{{sample}}_dge.summary.txt',species=config['META']['species'])
	output: 
		barcodes_species=expand('summary/{species}/{{sample}}_barcodes.csv', species=config['META']['species']),
		genes_pdf='plots/{sample}_species_plot_genes.pdf',
		transcripts_pdf='plots/{sample}_species_plot_transcripts.pdf'
	params:
		expected_cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells'])
	script: 
		'../scripts/plot_species_plot.R'