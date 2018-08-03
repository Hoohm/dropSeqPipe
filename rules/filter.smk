"""Filter data"""


#Which rules will be run on the host computer and not sent to nodes
localrules: count_reads, reads_after_trimming, plot_polyA_trim, plot_barcode_start_trim, plot_UMI_filtering, plot_CELL_filtering, plot_BC_drop, multiqc_trimmomatic

rule cutadapt_R1:
    input:
        R1='data/{sample}_R1.fastq.gz',
        adapters=config['FILTER']['cutadapt']['adapters-file']
    output:
        fastq=temp("data/{sample}/{sample}_trimmmed_R1.fastq.gz"),
        qc="logs/cutadapt/{sample}_R1.qc.txt"
    threads: 10
    conda: '../envs/cutadapt.yaml'
	shell:
		"""cutadapt -a file:{input.adapters} -a AAAAAAA -q 20 --discard-trimmed --minimum-length 20 --cores={threads} -o {output.fastq} {input.R1} > {output.qc}"""

rule cutadapt_R2:
    input:
        R1='data/{sample}_R2.fastq.gz',
        adapters=config['FILTER']['cutadapt']['adapters-file']
    output:
        fastq=temp("data/{sample}/{sample}_trimmmed_R2.fastq.gz"),
        qc="logs/cutadapt/{sample}_R2.qc.txt"
    threads: 10
    conda: '../envs/cutadapt.yaml'
	shell:
		"""cutadapt -a file:{input.adapters} -q 20 --discard-trimmed --minimum-length 15 --cores={threads} -o {output.fastq} {input.R1} > {output.qc}"""

rule repair:
	input:
		R1='data/{sample}/{sample}_trimmmed_R1.fastq.gz',
		R2='data/{sample}/{sample}_trimmmed_R2.fastq.gz'
	output:
		R1=temp('data/{sample}/{sample}_trimmmed_repaired_R1.fastq.gz'),
		R2=temp('data/{sample}/{sample}_trimmmed_repaired_R2.fastq.gz')
	conda: '../envs/bbmap.yaml'
	shell:
		"""repair.sh -Xmx400g in={input.R1} in2={input.R2} out1={output.R1} out2={output.R2} repair=t"""

rule fastq_to_sam:
	"""Create an empty bam file linking cell/UMI barcodes to reads"""
	input:
		R1='data/{sample}/{sample}_trimmmed_repaired_R1.fastq.gz',
		R2='data/{sample}/{sample}_trimmmed_repaired_R2.fastq.gz'
	output:
		temp('data/{sample}/{sample}_unaligned.bam')
	params:
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory'],
		picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar"
	conda: '../envs/picard.yaml'
	shell:
		"""java -Djava.io.tmpdir={params.temp_directory} -Xmx{params.memory} -jar {params.picard} FastqToSam\
		F1={input[0]}\
		F2={input[1]}\
		SM=DS O={output}"""


rule count_reads:
	input:
		'data/{sample}/{sample}_unaligned.bam'
	output:
		'logs/{sample}_reads_left.txt'
	conda: '../envs/samtools.yaml'
	shell:
		"""samtools view {input} | wc -l > {output}"""



rule sam_to_fastq:
	input:
		'data/{sample}/{sample}_unaligned.bam'
	params:
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory'],
		picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar"
	output:
		temp('data/{samples}/{sample}_trimmed_unmapped.fastq.gz')
	conda: '../envs/picard.yaml'
	shell:
		"""java -Xmx{params.memory} -jar -Djava.io.tmpdir={params.temp_directory}	{params.picard} SamToFastq\
		INPUT={input}\
		FASTQ=/dev/stdout COMPRESSION_LEVEL=0|\
		gzip > {output}"""

rule plot_BC_drop:
	input:
		Cell_tagged=expand('logs/{sample}_CELL_barcode.txt', sample=samples.index),
		UMI_tagged=expand('logs/{sample}_UMI_barcode.txt', sample=samples.index),
		reads_left=expand('logs/{sample}_reads_left.txt', sample=samples.index),
		trimmomatic_filtered=expand('logs/{sample}_reads_left_trim.txt', sample=samples.index)
	params:
		Cell_length=config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start']+1,
		UMI_length=config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start']+1,
		min_num_below_Cell=config['FILTER']['cell-barcode']['num-below-quality'],
		min_num_below_UMI=config['FILTER']['UMI-barcode']['num-below-quality'],
		min_Cell_quality=config['FILTER']['cell-barcode']['min-quality'],
		min_UMI_quality=config['FILTER']['UMI-barcode']['min-quality'],
		sample_names=lambda wildcards: samples.index,
		batches=lambda wildcards: samples.loc[samples.index, 'batch']
		
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/BC_drop.pdf'
	script:
		'../scripts/plot_BC_drop.R'

rule multiqc_trimmomatic:
	input:
		expand('logs/{sample}_trimlog.txt', sample=samples.index)
	params: '-m trimmomatic'
	output:
		html='reports/filter.html'
	wrapper:
		'0.21.0/bio/multiqc'