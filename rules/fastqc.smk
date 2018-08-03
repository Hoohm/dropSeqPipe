"""Get fastqc reports"""

#Which rules will be run on the host computer and not sent to nodes
localrules: multiqc_fastqc_reads, multiqc_fastqc_barcodes

rule fastqc_barcodes:
	"""Create fastqc report"""
	input: 
		'data/{sample}_R1.fastq.gz'
	output:
		html='logs/fastqc/{sample}_R1_fastqc.html',
		zip='logs/fastqc/{sample}_R1_fastqc.zip'
	params: '--extract'
	wrapper:
		'0.21.0/bio/fastqc'

rule fastqc_reads:
	"""Create fastqc report"""
	input: 
		'data/{sample}_R2.fastq.gz'
	output:
		html='logs/fastqc/{sample}_R2_fastqc.html',
		zip='logs/fastqc/{sample}_R2_fastqc.zip'
	params: '--extract'
	wrapper:
		'0.21.0/bio/fastqc'


rule multiqc_fastqc_barcodes:
	input:
		expand('logs/fastqc/{sample}_R1_fastqc.html', sample=samples.index)
	output:
		html='reports/fastqc_barcodes.html'
	params: '-m fastqc --ignore *_R2*'
	wrapper:
		'0.21.0/bio/multiqc'

rule multiqc_fastqc_reads:
	input: 
		expand('logs/fastqc/{sample}_R2_fastqc.html', sample=samples.index)
	output:
		html='reports/fastqc_reads.html'
	params: '-m fastqc --ignore *_R1*'
	wrapper:
		'0.21.0/bio/multiqc'