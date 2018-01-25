"""Extract expression fof single species"""

ruleorder: extract_umi_expression_whitelist > extract_umi_expression
ruleorder: extract_reads_expression_whitelist > extract_reads_expression
ruleorder: extract_umi_per_gene_whitelist > extract_umi_per_gene
ruleorder: SingleCellRnaSeqMetricsCollector_whitelist > SingleCellRnaSeqMetricsCollector
ruleorder: plot_rna_metrics_whitelist > plot_rna_metrics


rule extract_umi_expression:
	input:
		'data/{sample}_final.bam'
	output:
		'summary/{sample}_umi_expression_matrix.tsv'
	params:
		summary='summary/{sample}_dge.summary.txt',
		count_per_umi=config['EXTRACTION']['min_count_per_umi'],
		num_cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		bc_edit_distance=config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input}\
		O={output}\
		EDIT_DISTANCE={params.bc_edit_distance}\
		SUMMARY={params.summary}\
		NUM_CORE_BARCODES={params.num_cells}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""


rule extract_umi_expression_whitelist:
	input: 
		sample='data/{sample}_final.bam',
		barcode_whitelist='barcodes.csv'
	output:
		'summary/{sample}_umi_expression_matrix.tsv'
	params:
		summary='summary/{sample}_dge.summary.txt',
		count_per_umi=config['EXTRACTION']['min_count_per_umi'],
		bc_edit_distance=config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input.sample}\
		O={output}\
		EDIT_DISTANCE={params.bc_edit_distance}\
		CELL_BC_FILE={input.barcode_whitelist}\
		SUMMARY={params.summary}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""

rule extract_reads_expression_whitelist:
	input: 
		sample='data/{sample}_final.bam',
		barcode_whitelist='barcodes.csv'
	output:
		'summary/{sample}_counts_expression_matrix.tsv'
	params:
		summary='summary/{sample}_dge.summary.txt',
		count_per_umi=config['EXTRACTION']['min_count_per_umi'],
		bc_edit_distance=config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']		
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input.sample}\
		O={output}\
		EDIT_DISTANCE={params.bc_edit_distance}\
		OUTPUT_READS_INSTEAD=true\
		CELL_BC_FILE={input.barcode_whitelist}\
		SUMMARY={params.summary}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""

rule extract_reads_expression:
	input:
		'data/{sample}_final.bam'
	output:
		'summary/{sample}_counts_expression_matrix.tsv'
	params:
		count_per_umi=config['EXTRACTION']['min_count_per_umi'],
		num_cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		bc_edit_distance=config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']		
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DigitalExpression\
		I={input}\
		O={output}\
		EDIT_DISTANCE={params.bc_edit_distance}\
		OUTPUT_READS_INSTEAD=true\
		NUM_CORE_BARCODES={params.num_cells}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""


rule extract_umi_per_gene:
	input:
		'data/{sample}_final.bam'
	output:
		'logs/{sample}_umi_per_gene.tsv'
	params:
		num_cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		bc_edit_distance=config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']		
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p GatherMolecularBarcodeDistributionByGene\
		EDIT_DISTANCE={params.bc_edit_distance}\
		I={input}\
		O={output}\
		NUM_CORE_BARCODES={params.num_cells}"""

rule extract_umi_per_gene_whitelist:
	input: 
		sample='data/{sample}_final.bam',
		barcode_whitelist='barcodes.csv'
	output:
		'logs/{sample}_umi_per_gene.tsv'
	params:
		bc_edit_distance=config['EXTRACTION']['bc_edit_distance'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p GatherMolecularBarcodeDistributionByGene\
		EDIT_DISTANCE={params.bc_edit_distance}\
		I={input.sample}\
		O={output}\
		CELL_BC_FILE={input.barcode_whitelist}"""


rule SingleCellRnaSeqMetricsCollector:
	input: 
		refFlat=expand('{annotation_prefix}.refFlat', annotation_prefix=annotation_prefix),
		rRNA_intervals=expand('{reference_prefix}.rRNA.intervals', reference_prefix=reference_prefix),
		data='data/{sample}_final.bam'
	params:
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']
	output:
		'logs/{sample}_rna_metrics.txt'
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p SingleCellRnaSeqMetricsCollector\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		NUM_CORE_BARCODES={params.cells}\
		RIBOSOMAL_INTERVALS={input.rRNA_intervals}
		"""

rule SingleCellRnaSeqMetricsCollector_whitelist:
	input:
		sample='data/{sample}_final.bam',
		barcode_whitelist='barcodes.csv'
	params:
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		refFlat=expand('{annotation_prefix}.refFlat', annotation_prefix=annotation_prefix),
		rRNA_intervals=expand('{reference_prefix}.rRNA.intervals', reference_prefix=reference_prefix),
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		TMPDIR=config['LOCAL']['TMPDIR'],
		memory=config['LOCAL']['MEMORY']
	output:
		'logs/{sample}_rna_metrics.txt'
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p SingleCellRnaSeqMetricsCollector\
		INPUT={input.sample}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={params.refFlat}\
		CELL_BC_FILE={input.barcode_whitelist}\
		RIBOSOMAL_INTERVALS={params.rRNA_intervals}
		"""

rule plot_umi_per_gene:
	input:
		expand('logs/{sample}_umi_per_gene.tsv',sample=samples.index)
	output:
		pdf='plots/umi_per_gene_distribution.pdf',
		png='plots/png/umi_per_gene_distribution.png'
	params:
		sample_names=lambda wildcards: samples.index	
	script:
		"../scripts/plot_umi_distribution.R"

rule plot_rna_metrics:
	input:
		'logs/{sample}_rna_metrics.txt'
	output:
		pdf='plots/{sample}_rna_metrics.pdf',
		png='plots/png/{sample}_rna_metrics.png'
	params: 
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells']
	script:
		'../scripts/plot_rna_metrics.R'

rule plot_rna_metrics_whitelist:
	input:
		rna_metrics='logs/{sample}_rna_metrics.txt',
		barcodes='barcodes.csv'
	output:
		pdf='plots/{sample}_rna_metrics.pdf',
		png='plots/png/{sample}_rna_metrics.png'
	script:
		'../scripts/plot_rna_metrics.R'

rule merge_umi:
	input:
		expand('summary/{sample}_umi_expression_matrix.tsv', sample=samples.index)
	params:
		sample_names=lambda wildcards: samples.index
	output:
		'summary/umi_expression_matrix.tsv'
	script:
		"../scripts/merge_counts_single.R"

rule merge_counts:
	input:
		expand('summary/{sample}_counts_expression_matrix.tsv', sample=samples.index)
	params:
		sample_names=lambda wildcards: samples.index
	output:
		'summary/counts_expression_matrix.tsv'
	script:
		"../scripts/merge_counts_single.R"
