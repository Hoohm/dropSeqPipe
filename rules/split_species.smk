"""Extract species specific expression to prepare the species plot."""

ruleorder: extract_all_umi_expression_whitelist_species > extract_all_umi_expression_species


rule split_bam_species:
	input: 'data/{sample}_final.bam'
	output: 'data/{sample}_{species}_unfiltered.bam'
	params:
		species= lambda wildcards: wildcards.species,
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory=config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p FilterBAM\
		REF_SOFT_MATCHED_RETAINED={params.species}\
		INPUT={input}\
		OUTPUT={output}"""


rule extract_all_umi_expression_species:
	input: 'data/{sample}_{species}_unfiltered.bam'
	output:
		umi_matrix = temp('summary/{sample}_{species}_unfiltered_umi_expression_matrix.tsv'),
		summary = 'summary/{sample}_{species}_dge.summary.txt'

	params:
		count_per_umi = config['EXTRACTION']['min_count_per_umi'],
		num_cells = lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		bc_edit_distance = config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory=config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input}\
		O={output.umi_matrix}\
		SUMMARY={output.summary}\
		EDIT_DISTANCE={params.bc_edit_distance}\
		NUM_CORE_BARCODES={params.num_cells}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""

rule extract_all_umi_expression_whitelist_species:
	input: 
		sample = 'data/{sample}_{species}_unfiltered.bam',
		barcode_whitelist = 'barcodes.csv'
	output:
		umi_matrix = temp('summary/{sample}_{species}_unfiltered_umi_expression_matrix.tsv'),
		summary = 'summary/{sample}_{species}_dge.summary.txt'
	params:
		count_per_umi = config['EXTRACTION']['min_count_per_umi'],
		bc_edit_distance = config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory=config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input.sample}\
		O={output.umi_matrix}\
		SUMMARY={output.summary}\
		EDIT_DISTANCE={params.bc_edit_distance}\
		CELL_BC_FILE={input.barcode_whitelist}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""


rule plot_barnyard:
	input:
		expand('summary/{{sample}}_{species}_dge.summary.txt',species=config['META']['species'])
	output: 
		barcodes_species = expand('summary/{{sample}}_{species}_barcodes.csv', species=config['META']['species']),
		barnyard_plot_genes = 'plots/{sample}_species_plot_genes.pdf',
		barnyard_plot_transcripts = 'plots/{sample}_species_plot_transcripts.pdf'
	params:
		expected_cells = lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells'])
	script: 
		'../scripts/plot_species_plot.R'

# rule species_plot:
# 	input: expand('summary/{{sample}}_{species}_dge.summary.txt',species=config['META']['species'])
# 	output:
# 		plot = 'plots/{sample}_species_plot.pdf'
# 	params:
# 		sample_names = lambda wildcards: samples.index
# 	script: 
# 		'../scripts/plot_species_plot.R'