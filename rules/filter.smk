"""Filter data"""


#Which rules will be run on the host computer and not sent to nodes
localrules: count_reads, reads_after_trimming, plot_polyA_trim, plot_barcode_start_trim, plot_UMI_filtering, plot_CELL_filtering, plot_BC_drop, multiqc_trimmomatic

rule fastq_to_sam:
	"""Create an empty bam file linking cell/UMI barcodes to reads"""
	input:
		R1='data/{sample}_R1.fastq.gz',
		R2='data/{sample}_R2.fastq.gz'
	output:
		temp('data/{sample}_unaligned.bam')
	params:
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory'],
		picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar"
	conda: '../envs/picard.yaml'
	shell:
		"""java -Djava.io.tmpdir={params.temp_directory} -Xmx{params.memory} -jar {params.picard} FastqToSam\
		F1={input.R1}\
		F2={input.R2}\
		SM=DS O={output}"""


rule BC_tags:
	input:
		'data/{sample}_unaligned.bam'
	output: 
		data=temp('data/{sample}_BC_tagged_unmapped.bam'),
		BC_summary='logs/{sample}_CELL_barcode.txt'
	params:
		BC_start=config['FILTER']['cell-barcode']['start'],
		BC_end=config['FILTER']['cell-barcode']['end'],
		BC_minQuality=config['FILTER']['cell-barcode']['min-quality'],
		BC_minQuality_num=config['FILTER']['cell-barcode']['num-below-quality']+1,
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && TagBamWithReadSequenceExtended -m {params.memory}\
		SUMMARY={output.BC_summary}\
		BASE_RANGE={params.BC_start}-{params.BC_end}\
		BASE_QUALITY={params.BC_minQuality}\
		BARCODED_READ=1\
		DISCARD_READ=false\
		TAG_NAME=XC\
		NUM_BASES_BELOW_QUALITY={params.BC_minQuality_num}\
		INPUT={input}\
		OUTPUT={output.data}
		"""

rule UMI_tags:
	input:
		'data/{sample}_BC_tagged_unmapped.bam'
	output:
		bam=temp('data/{sample}_BC_UMI_tagged_unmapped.bam'),
		UMI_summary='logs/{sample}_UMI_barcode.txt'
	params:
		UMI_start=config['FILTER']['UMI-barcode']['start'],
		UMI_end=config['FILTER']['UMI-barcode']['end'],
		UMI_minQuality=config['FILTER']['UMI-barcode']['min-quality'],
		UMI_minQuality_num=config['FILTER']['UMI-barcode']['num-below-quality']+1,
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && TagBamWithReadSequenceExtended -m {params.memory}\
		SUMMARY={output.UMI_summary}\
		BASE_RANGE={params.UMI_start}-{params.UMI_end}\
		BASE_QUALITY={params.UMI_minQuality}\
		BARCODED_READ=1\
		DISCARD_READ=true\
		TAG_NAME=XM\
		NUM_BASES_BELOW_QUALITY={params.UMI_minQuality_num}\
		INPUT={input}\
		OUTPUT={output.bam}
		"""
rule filter_tags:
	input:
		'data/{sample}_BC_UMI_tagged_unmapped.bam'
	output: 
		temp('data/{sample}_tags_filtered_unmapped.bam')
	params:		
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} &&  FilterBAM -m {params.memory}\
		TAG_REJECT=XQ\
		INPUT={input}\
		OUTPUT={output}
		"""

rule count_reads:
	input:
		'data/{sample}_tags_filtered_unmapped.bam'
	output:
		'logs/{sample}_reads_left.txt'
	conda: '../envs/samtools.yaml'
	shell:
		"""samtools view {input} | wc -l > {output}"""


rule start_trim:
	input:
		'data/{sample}_tags_filtered_unmapped.bam'
	output:
		data=temp('data/{sample}_tags_start_filtered_unmapped.bam'),
		trim_summary='logs/{sample}_start_trim.txt'
	params:
		SmartAdapter=config['FILTER']['5-prime-smart-adapter'],
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && TrimStartingSequence -m {params.memory}\
		OUTPUT_SUMMARY={output.trim_summary}\
		SEQUENCE={params.SmartAdapter}\
		MISMATCHES=1\
		NUM_BASES=6\
		INPUT={input}\
		OUTPUT={output.data}
		"""
rule polya_trim:
	input:
		'data/{sample}_tags_start_filtered_unmapped.bam'
	output:
		data='data/{sample}_trimmed_unmapped.bam',
		trim_summary='logs/{sample}_polyA_trim.txt'
	params:		
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && PolyATrimmer -m {params.memory}\
		OUTPUT_SUMMARY={output.trim_summary}\
		MISMATCHES=0\
		NUM_BASES=5\
		INPUT={input}\
		OUTPUT={output.data}
		"""

rule sam_to_fastq:
	input:
		'data/{sample}_trimmed_unmapped.bam'
	params:
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory'],
		picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar"
	output:
		temp('data/{sample}_trimmed_unmapped.fastq.gz')
	conda: '../envs/picard.yaml'
	shell:
		"""java -Xmx{params.memory} -jar -Djava.io.tmpdir={params.temp_directory}	{params.picard} SamToFastq\
		INPUT={input}\
		FASTQ=/dev/stdout COMPRESSION_LEVEL=0|\
		gzip > {output}"""



rule trim_single:
	input:
		'data/{sample}_trimmed_unmapped.fastq.gz'
	output:
		data='data/{sample}_filtered.fastq.gz'
	log:
		'logs/{sample}_trimlog.txt'
	params:
		trimmer=['LEADING:{} TRAILING:{} SLIDINGWINDOW:{}:{} MINLEN:{} ILLUMINACLIP:{}:{}:{}:{}'.format(
			config['FILTER']['trimmomatic']['LEADING'],
			config['FILTER']['trimmomatic']['TRAILING'],
			config['FILTER']['trimmomatic']['SLIDINGWINDOW']['windowSize'],
			config['FILTER']['trimmomatic']['SLIDINGWINDOW']['requiredQuality'],
			config['FILTER']['trimmomatic']['MINLEN'],
			config['FILTER']['trimmomatic']['adapters-file'],
			config['FILTER']['trimmomatic']['ILLUMINACLIP']['seedMismatches'],
			config['FILTER']['trimmomatic']['ILLUMINACLIP']['palindromeClipThreshold'],
			config['FILTER']['trimmomatic']['ILLUMINACLIP']['simpleClipThreshold'])],
		extra='-threads 10'
	threads: 10
	wrapper:
		'0.27.1/bio/trimmomatic/se'


rule reads_after_trimming:
	input:
		'data/{sample}_filtered.fastq.gz'
	output:
		'logs/{sample}_reads_left_trim.txt'
	conda: '../envs/samtools.yaml'
	shell:
		"""zcat {input} | wc -l > {output}"""

rule plot_polyA_trim:
	input:
		'logs/{sample}_polyA_trim.txt'
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/{sample}_polya_trimmed.pdf'
	script:
		'../scripts/plot_polyA_trim.R'

rule plot_barcode_start_trim:
	input:
		'logs/{sample}_start_trim.txt'
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/{sample}_start_trim.pdf'
	script:
		'../scripts/plot_start_trim.R'


rule plot_UMI_filtering:
	input:
		'logs/{sample}_UMI_barcode.txt'
	params: 
		min_quality=config['FILTER']['UMI-barcode']['min-quality'],
		num_below_quality=config['FILTER']['UMI-barcode']['num-below-quality']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/{sample}_UMI_dropped.pdf'
	script:
		'../scripts/plot_umi_drop.R'

rule plot_CELL_filtering:
	input:
		'logs/{sample}_CELL_barcode.txt'
	params:
		min_quality=config['FILTER']['cell-barcode']['min-quality'],
		num_below_quality=config['FILTER']['cell-barcode']['num-below-quality']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/{sample}_CELL_dropped.pdf'

	script:
		'../scripts/plot_cell_drop.R'

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
		'0.27.1/bio/multiqc'