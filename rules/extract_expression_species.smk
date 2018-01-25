"""Extract expression fof mixed species"""

rule extract_umi_expression_species:
	input:
		data='data/{sample}_{species}_unfiltered.bam',
		barcode_whitelist='summary/{sample}_{species}_barcodes.csv'
	output:
		'summary/{sample}_{species}_umi_expression_matrix.txt'
	params:
		count_per_umi=config['EXTRACTION']['min_count_per_umi'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory=config['LOCAL']['MEMORY'],
		TMPDIR=config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input.data}\
		O={output}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}\
		CELL_BC_FILE={input.barcode_whitelist}"""

rule extract_reads_expression_species:
	input:
		data='data/{sample}_{species}_unfiltered.bam',
		barcode_whitelist='summary/{sample}_{species}_barcodes.csv'
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory=config['LOCAL']['MEMORY'],
		TMPDIR=config['LOCAL']['TMPDIR']
	output:
		'summary/{sample}_{species}_counts_expression_matrix.txt'
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input.data}\
		O={output}\
		CELL_BC_FILE={input.barcode_whitelist}\
		OUTPUT_READS_INSTEAD=true"""


rule extract_umi_per_gene_species:
	input:
		data='data/{sample}_{species}_unfiltered.bam',
		barcode_whitelist='summary/{sample}_{species}_barcodes.csv'
	params:	
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory=config['LOCAL']['MEMORY'],
		TMPDIR=config['LOCAL']['TMPDIR']
	output:
		'logs/{sample}_{species}_umi_per_gene.tsv'
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p GatherMolecularBarcodeDistributionByGene\
		I={input.data}\
		O={output}\
		CELL_BC_FILE={input.barcode_whitelist}"""

rule SingleCellRnaSeqMetricsCollector_whitelist_species:
	input:
		sample='data/{sample}_{species}_unfiltered.bam',
		barcode_whitelist='summary/{sample}_{species}_barcodes.csv'
	params:
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		refFlat=expand('{annotation_prefix}.refFlat', annotation_prefix=annotation_prefix),
		rRNA_intervals=expand('{reference_prefix}.rRNA.intervals', reference_prefix=reference_prefix),
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory=config['LOCAL']['MEMORY'],
		TMPDIR=config['LOCAL']['TMPDIR']
	output:
		'logs/{sample}_{species}_rna_metrics.txt'
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p SingleCellRnaSeqMetricsCollector\
		INPUT={input.sample}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={params.refFlat}\
		CELL_BC_FILE={input.barcode_whitelist}\
		RIBOSOMAL_INTERVALS={params.rRNA_intervals}
		"""
rule plot_rna_metrics_species:
	input:
		rna_metrics='logs/{sample}_{species}_rna_metrics.txt',
		barcode='summary/{sample}_{species}_barcodes.csv'
	output:
		pdf='plots/species_{sample}_{species}_rna_metrics.pdf',
		png='plots/png/species_{sample}_{species}_rna_metrics.png'
	script:
		'../scripts/plot_rna_metrics.R'


rule merge_umi_species:
	input:
		expand('summary/{sample}_{{species}}_umi_expression_matrix.txt', sample=samples.index)
	params:
		sample_names=lambda wildcards: samples.index
	output:
		'summary/Experiment_{species}_umi_expression_matrix.tsv'
	script:
		"../scripts/merge_counts.R"

rule merge_counts_species:
	input:
		expand('summary/{sample}_{{species}}_counts_expression_matrix.txt', sample=samples.index)
	params:
		sample_names=lambda wildcards: samples.index
	output:
		'summary/Experiment_{species}_counts_expression_matrix.tsv'
	script:
		"../scripts/merge_counts.R"