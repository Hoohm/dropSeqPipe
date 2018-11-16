

ruleorder: extend_barcode_whitelist > get_top_barcodes
localrules: extend_barcode_whitelist, extend_barcode_top

rule get_top_barcodes:
	input:
		"data/{sample}/trimmmed_repaired_R1.fastq.gz"
	output:
		'data/{sample}/top_barcodes.csv'
	conda: '../envs/umi_tools.yaml'
	params:
		cell_barcode_length=(config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1),
		umi_barcode_length=(config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'] + 1),
		num_cells=lambda wildcards: round(int(samples.loc[wildcards.sample,'expected_cells'])*1.2),
	shell:
		"""umi_tools whitelist\
		--stdin {input}\
		--bc-pattern='(?P<cell_1>.{{{params.cell_barcode_length}}})(?P<umi_1>.{{{params.umi_barcode_length}}})'\
		--extract-method=regex\
		--set-cell-number={params.num_cells}\
		--log2stderr > {output}"""

rule get_cell_whitlist:
	input:
		'data/{sample}/top_barcodes.csv'
	output:
		'data/{sample}/barcodes.csv'
	shell:
		"""cat {input} | cut -f 1 > {output}"""

rule extend_barcode_whitelist:
	input:
		whitelist='barcodes.csv'
	output:
		barcodes='data/{sample}/barcodes.csv',
		barcode_ref='data/{sample}/barcode_ref.pkl',
		barcode_ext_ref='data/{sample}/barcode_ext_ref.pkl',
		barcode_mapping='data/{sample}/empty_barcode_mapping.pkl'
	script:
		'../scripts/generate_extended_ref.py'



rule extend_barcode_top:
	input:
		whitelist='data/{sample}/top_barcodes.csv'
	output:
		barcode_ref='data/{sample}/barcode_ref.pkl',
		barcode_ext_ref='data/{sample}/barcode_ext_ref.pkl',
		barcode_mapping='data/{sample}/empty_barcode_mapping.pkl'
	script:
		'../scripts/umi_tools_extended_ref.py'


rule repair_barcodes:
	input:
		bam='data/{sample}/Aligned.merged.bam',
		barcode_ref='data/{sample}/barcode_ref.pkl',
		barcode_ext_ref='data/{sample}/barcode_ext_ref.pkl',
		barcode_mapping='data/{sample}/empty_barcode_mapping.pkl'
	conda: '../envs/merge_bam.yaml'
	output:
		bam=temp('data/{sample}/Aligned.repaired.bam'),
		barcode_mapping_counts='data/{sample}/barcode_mapping_counts.pkl'
	script:
		'../scripts/repair_barcodes.py'