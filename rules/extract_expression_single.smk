"""Extract expression fof single species"""

ruleorder: extract_umi_expression_whitelist > extract_umi_expression
ruleorder: extract_reads_expression_whitelist > extract_reads_expression
ruleorder: extract_umi_per_gene_whitelist > extract_umi_per_gene
ruleorder: SingleCellRnaSeqMetricsCollector_whitelist > SingleCellRnaSeqMetricsCollector
ruleorder: plot_rna_metrics_whitelist > plot_rna_metrics

#Which rules will be run on the host computer and not sent to nodes
localrules: plot_umi_per_gene, plot_rna_metrics, plot_rna_metrics_whitelist, merge_umi, merge_counts, violine_plots

rule extract_umi_expression:
	input:
		'data/{sample}/final.bam'
	output:
		'data/{sample}/umi_expression_matrix.tsv'
	log:
		summary='data/{sample}/dge.summary.txt',
	params:
		count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
		num_cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
		I={input}\
		O={output}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		SUMMARY={log.summary}\
		NUM_CORE_BARCODES={params.num_cells}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""


rule extract_umi_expression_whitelist:
	input: 
		data='data/{sample}/final.bam',
		barcode_whitelist='barcodes.csv'
	output:
		'data/{sample}/umi_expression_matrix.tsv'
	log:
		summary='data/{sample}/dge.summary.txt',
	params:
		count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],		
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
		I={input.data}\
		O={output}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		CELL_BC_FILE={input.barcode_whitelist}\
		SUMMARY={log.summary}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""

rule extract_reads_expression_whitelist:
	input: 
		data='data/{sample}/final.bam',
		barcode_whitelist='barcodes.csv'
	output:
		'data/{sample}/counts_expression_matrix.tsv'
	params:
		summary='data/{sample}/dge.summary.txt',
		count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],		
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	conda: '../envs/dropseq_tools.yaml'	
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
		I={input.data}\
		O={output}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		OUTPUT_READS_INSTEAD=true\
		CELL_BC_FILE={input.barcode_whitelist}\
		SUMMARY={params.summary}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""

rule extract_reads_expression:
	input:
		'data/{sample}/final.bam'
	output:
		'data/{sample}/counts_expression_matrix.tsv'
	params:
		count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
		num_cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],	
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	conda: '../envs/dropseq_tools.yaml'	
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DigitalExpression -m {params.memory}\
		I={input}\
		O={output}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		OUTPUT_READS_INSTEAD=true\
		NUM_CORE_BARCODES={params.num_cells}\
		MIN_BC_READ_THRESHOLD={params.count_per_umi}"""


rule extract_umi_per_gene:
	input:
		'data/{sample}/final.bam'
	output:
		'logs/dropseq_tools/{sample}_umi_per_gene.tsv'
	params:
		num_cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],	
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	conda: '../envs/dropseq_tools.yaml'	
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && GatherMolecularBarcodeDistributionByGene -m {params.memory}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		I={input}\
		O={output}\
		NUM_CORE_BARCODES={params.num_cells}"""

rule extract_umi_per_gene_whitelist:
	input: 
		data='data/{sample}/final.bam',
		barcode_whitelist='barcodes.csv'
	output:
		'logs/dropseq_tools/{sample}_umi_per_gene.tsv'
	params:
		cellBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],	
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && GatherMolecularBarcodeDistributionByGene -m {params.memory}\
		EDIT_DISTANCE={params.cellBarcodeEditDistance}\
		I={input.data}\
		O={output}\
		CELL_BC_FILE={input.barcode_whitelist}"""


rule SingleCellRnaSeqMetricsCollector:
	input: 
		data='data/{sample}/final.bam',
		refFlat="{}.refFlat".format(annotation_prefix),
		rRNA_intervals="{}.rRNA.intervals".format(reference_prefix),
	params:
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],		
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	output:
		'logs/dropseq_tools/{sample}_rna_metrics.txt'
	wildcard_constraints:
		sample="({})".format("|".join(samples.index))
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && SingleCellRnaSeqMetricsCollector -m {params.memory}\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		NUM_CORE_BARCODES={params.cells}\
		RIBOSOMAL_INTERVALS={input.rRNA_intervals}
		"""

rule SingleCellRnaSeqMetricsCollector_whitelist:
	input:
		data='data/{sample}/final.bam',
		barcode_whitelist='barcodes.csv',
		refFlat="{}.refFlat".format(annotation_prefix),
		rRNA_intervals="{}.rRNA.intervals".format(reference_prefix)
	params:		
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	output:
		'logs/dropseq_tools/{sample}_rna_metrics.txt'
	wildcard_constraints:
		sample="({})".format("|".join(samples.index))
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && SingleCellRnaSeqMetricsCollector -m {params.memory}\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		CELL_BC_FILE={input.barcode_whitelist}\
		RIBOSOMAL_INTERVALS={input.rRNA_intervals}
		"""

rule plot_umi_per_gene:
	input:
		expand('logs/dropseq_tools/{sample}_umi_per_gene.tsv',sample=samples.index)
	params:
		sample_names=lambda wildcards: samples.index
	conda: '../envs/plots.yaml'	
	output:
		pdf='plots/umi_per_gene_distribution.pdf'
	script:
		"../scripts/plot_umi_distribution.R"

rule plot_rna_metrics:
	input:
		'logs/dropseq_tools/{sample}_rna_metrics.txt'
	params: 
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/rna_metrics/{sample}_rna_metrics.pdf'
	script:
		'../scripts/plot_rna_metrics.R'

rule plot_rna_metrics_whitelist:
	input:
		rna_metrics='logs/dropseq_tools/{sample}_rna_metrics.txt',
		barcodes='barcodes.csv'
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/rna_metrics/{sample}_rna_metrics.pdf',
	wildcard_constraints:
		sample="({})".format("|".join(samples.index))
	script:
		'../scripts/plot_rna_metrics.R'

rule merge_umi:
	input:
		expand('data/{sample}/umi_expression_matrix.tsv', sample=samples.index)
	params:
		sample_names=lambda wildcards: samples.index
	conda: '../envs/merge.yaml'
	output:
		'summary/umi_expression_matrix.tsv'
	script:
		"../scripts/merge_counts_single.R"

rule merge_counts:
	input:
		expand('data/{sample}/counts_expression_matrix.tsv', sample=samples.index)
	params:
		sample_names=lambda wildcards: samples.index
	conda: '../envs/merge.yaml'
	output:
		'summary/counts_expression_matrix.tsv'
	script:
		"../scripts/merge_counts_single.R"

rule violine_plots:
	input:
		UMIs='summary/umi_expression_matrix.tsv',
		counts='summary/counts_expression_matrix.tsv',
		design='samples.csv'
	conda: '../envs/plots_ext.yaml'
	output:
		pdf_violine='plots/violinplots_comparison_UMI.pdf',
		pdf_umivscounts='plots/UMI_vs_counts.pdf',
		pdf_umi_vs_gene='plots/UMI_vs_gene.pdf',
		pdf_count_vs_gene='plots/Count_vs_gene.pdf',
		R_objects='summary/R_Seurat_objects.rdata'
	script:
		'../scripts/plot_violine.R'