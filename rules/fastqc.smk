"""Get fastqc reports"""

def get_fastq(wildcards):
    return 'data/' + samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

rule fastqc:
	"""Create fastqc report"""
	input: 
		R1 = 'data/{sample}_R1.fastq.gz',
		R2 = 'data/{sample}_R2.fastq.gz'
	output: 'logs/{sample}_R1_fastqc.html'
	conda: '../envs/fastqc.yaml'
	threads: 2
	shell:
		"""fastqc {input.R1} {input.R2} -t 2 -o logs --extract"""

rule multiqc_fastqc:
	input: expand('logs/{sample}_R1_fastqc.html', sample=samples.index)
	output:
		html = 'reports/fastqc.html',
		txt = 'reports/fastqc_data/multiqc_general_stats.txt'
	shell:
		"""multiqc logs/ -m fastqc -o reports/ -n fastqc.html -f"""