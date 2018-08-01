"""Extract expression fof mixed species"""

#Which rules will be run on the host computer and not sent to nodes
localrules: plot_rna_metrics_species, merge_umi_species, merge_counts_species

rule extract_umi_expression_species:
	input:
		data='data/{species}/{sample}_unfiltered.bam',
		barcode_whitelist='summary/{species}/{sample}_barcodes.csv'
	output:
		'summary/{species}/{sample}_umi_expression_matrix.txt'
	params:
		count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""{params.dropseq_wrapper} -t {params.temp_directory} -m {params.memory} -p DigitalExpression\
		I={input.data}\
		O={output}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}\
		CELL_BC_FILE={input.barcode_whitelist}"""

rule extract_reads_expression_species:
	input:
		data='data/{species}/{sample}_unfiltered.bam',
		barcode_whitelist='summary/{species}/{sample}_barcodes.csv'
	params:
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'summary/{species}/{sample}_counts_expression_matrix.txt'
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""{params.dropseq_wrapper} -t {params.temp_directory} -m {params.memory} -p DigitalExpression\
		I={input.data}\
		O={output}\
		CELL_BC_FILE={input.barcode_whitelist}\
		OUTPUT_READS_INSTEAD=true"""


rule extract_umi_per_gene_species:
	input:
		data='data/{species}/{sample}_unfiltered.bam',
		barcode_whitelist='summary/{species}/{sample}_barcodes.csv'
	params:	
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'logs/{species}/{sample}_umi_per_gene.tsv'
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""{params.dropseq_wrapper} -t {params.temp_directory} -m {params.memory} -p GatherMolecularBarcodeDistributionByGene\
		I={input.data}\
		O={output}\
		CELL_BC_FILE={input.barcode_whitelist}"""

rule SingleCellRnaSeqMetricsCollector_whitelist_species:
	input:
		data='data/{species}/{sample}_unfiltered.bam',
		barcode_whitelist='summary/{species}/{sample}_barcodes.csv',
		refFlat='{}.refFlat'.format(annotation_prefix),
		rRNA_intervals='{}.rRNA.intervals'.format(reference_prefix)
	params:
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'logs/{species}/{sample}_rna_metrics.txt'
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""{params.dropseq_wrapper} -t {params.temp_directory} -m {params.memory} -p SingleCellRnaSeqMetricsCollector\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		CELL_BC_FILE={input.barcode_whitelist}\
		RIBOSOMAL_INTERVALS={input.rRNA_intervals}
		"""
rule plot_rna_metrics_species:
	input:
		rna_metrics='logs/{species}/{sample}_rna_metrics.txt',
		barcode='summary/{species}/{sample}_barcodes.csv'
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/{species}/{sample}_rna_metrics.pdf'
	script:
		'../scripts/plot_rna_metrics.R'


rule merge_umi_species:
	input:
		expand('summary/{{species}}/{sample}_umi_expression_matrix.txt', sample=samples.index)
	conda: '../envs/merge.yaml'
	output:
		'summary/Experiment_{species}_umi_expression_matrix.tsv'
	params:
		sample_names=lambda wildcards: samples.index
	script:
		"../scripts/merge_counts.R"

rule merge_counts_species:
	input:
		expand('summary/{{species}}/{sample}_counts_expression_matrix.txt', sample=samples.index)
	conda: '../envs/merge.yaml'
	output:
		'summary/Experiment_{species}_counts_expression_matrix.tsv'
	params:
		sample_names=lambda wildcards: samples.index
	script:
		"../scripts/merge_counts.R"