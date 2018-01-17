"""Filter data"""

rule fastq_to_sam:
	"""Create an empty bam file linking cell/UMI barcodes to reads"""
	input:
		R1 = 'data/{sample}_R1.fastq.gz',
		R2 = 'data/{sample}_R2.fastq.gz'
	output:
		temp('data/{sample}_unaligned.bam')
	params:
		TMPDIR = config['LOCAL']['TMPDIR'],
		memory = config['LOCAL']['MEMORY'],
		picard = "$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar"
	shell:
		"""java -Djava.io.tmpdir={params.TMPDIR} -Xmx{params.memory} -jar {params.picard} FastqToSam\
		F1={input[0]}\
		F2={input[1]}\
		SM=DS O={output}"""

rule BC_tags:
	input:
		'data/{sample}_unaligned.bam'
	output: 
		bam = temp('data/{sample}_BC_tagged_unmapped.bam'),
		BC_summary = 'logs/{sample}_CELL_barcode.txt'
	params:
		BC_start = config['FILTER']['Cell_barcode']['start'],
		BC_end = config['FILTER']['Cell_barcode']['end'],
		BC_min_quality = config['FILTER']['Cell_barcode']['min_quality'],
		BC_min_quality_num = config['FILTER']['Cell_barcode']['num_below_quality']+1,
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p TagBamWithReadSequenceExtended\
		SUMMARY={output.BC_summary}\
		BASE_RANGE={params.BC_start}-{params.BC_end}\
		BASE_QUALITY={params.BC_min_quality}\
		BARCODED_READ=1\
		DISCARD_READ=false\
		TAG_NAME=XC\
		NUM_BASES_BELOW_QUALITY={params.BC_min_quality_num}\
		INPUT={input}\
		OUTPUT={output.bam}
		"""

rule UMI_tags:
	input: 'data/{sample}_BC_tagged_unmapped.bam'
	output:
		bam = temp('data/{sample}_BC_UMI_tagged_unmapped.bam'),
		UMI_summary = 'logs/{sample}_UMI_barcode.txt'
	params:
		UMI_start = config['FILTER']['UMI']['start'],
		UMI_end = config['FILTER']['UMI']['end'],
		UMI_min_quality = config['FILTER']['UMI']['min_quality'],
		UMI_min_quality_num = config['FILTER']['UMI']['num_below_quality']+1,
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p TagBamWithReadSequenceExtended\
		SUMMARY={output.UMI_summary}\
		BASE_RANGE={params.UMI_start}-{params.UMI_end}\
		BASE_QUALITY={params.UMI_min_quality}\
		BARCODED_READ=1\
		DISCARD_READ=true\
		TAG_NAME=XM\
		NUM_BASES_BELOW_QUALITY={params.UMI_min_quality_num}\
		INPUT={input}\
		OUTPUT={output.bam}
		"""
rule filter_tags:
	input: 'data/{sample}_BC_UMI_tagged_unmapped.bam'
	output: 
		temp('data/{sample}_tags_filtered_unmapped.bam')
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p FilterBAM\
		TAG_REJECT=XQ\
		INPUT={input}\
		OUTPUT={output}
		"""

rule count_reads:
	input: 'data/{sample}_tags_filtered_unmapped.bam'
	output: 'logs/{sample}_reads_left.txt'
	shell:
		"""samtools view {input} | wc -l > {output}"""


rule start_trim:
	input: 'data/{sample}_tags_filtered_unmapped.bam'
	output:
		bam = temp('data/{sample}_tags_start_filtered_unmapped.bam'),
		trim_summary = 'logs/{sample}_start_trim.txt'
	params:
		SmartAdapter = config['FILTER']['5PrimeSmartAdapter'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p TrimStartingSequence\
		OUTPUT_SUMMARY={output.trim_summary}\
		SEQUENCE={params.SmartAdapter}\
		MISMATCHES=1\
		NUM_BASES=6\
		INPUT={input}\
		OUTPUT={output.bam}
		"""
rule polya_trim:
	input: 'data/{sample}_tags_start_filtered_unmapped.bam'
	output:
		bam ='data/{sample}_trimmed_unmapped.bam',
		trim_summary = 'logs/{sample}_polyA_trim.txt'
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p PolyATrimmer\
		OUTPUT_SUMMARY={output.trim_summary}\
		MISMATCHES=0\
		NUM_BASES=5\
		INPUT={input}\
		OUTPUT={output.bam}
		"""

rule sam_to_fastq:
	input: 'data/{sample}_trimmed_unmapped.bam'
	output: temp('data/{sample}_trimmed_unmapped.fastq.gz')
	params:
		TMPDIR = config['LOCAL']['TMPDIR'],
		memory = config['LOCAL']['MEMORY'],
		picard = "$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar"
	shell:
		"""java -Xmx{params.memory} -jar -Djava.io.tmpdir={params.TMPDIR}	{params.picard} SamToFastq\
		INPUT={input}\
		FASTQ=/dev/stdout COMPRESSION_LEVEL=0|\
		gzip > {output}"""

rule trim_single:
	input: 'data/{sample}_trimmed_unmapped.fastq.gz'
	output:
		data = 'data/{sample}_filtered.fastq.gz',
		log = 'logs/{sample}_trimlog.txt'
	params:
		ILLUMINACLIP = config['FILTER']['IlluminaClip']
	threads: 2
	shell:
		"""trimmomatic\
		SE {input} {output.data}\
		-threads {threads}\
		ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/{params.ILLUMINACLIP}:2:30:10\
		LEADING:3\
		TRAILING:3\
		SLIDINGWINDOW:4:20\
		MINLEN:15 > {output.log} 2>&1"""

rule plot_polyA_trim:
	input: 'logs/{sample}_polyA_trim.txt'
	output:
		pdf = 'plots/{sample}_polya_trimmed.pdf',
		png = 'plots/png/{sample}_polya_trimmed.png'
	script:
		'../scripts/plot_polyA_trim.R'

rule plot_barcode_start_trim:
	input: 'logs/{sample}_start_trim.txt'
	output:
		pdf = 'plots/{sample}_start_trim.pdf',
		png = 'plots/png/{sample}_start_trim.png'
	script:
		'../scripts/plot_start_trim.R'


rule plot_UMI_filtering:
	input: 'logs/{sample}_UMI_barcode.txt'
	output:
		pdf = 'plots/{sample}_UMI_dropped.pdf',
		png = 'plots/png/{sample}_UMI_dropped.png'
	params: 
		min_quality = config['FILTER']['UMI']['min_quality'],
		num_below_quality = config['FILTER']['UMI']['num_below_quality']
	script:
		'../scripts/plot_umi_drop.R'

rule plot_CELL_filtering:
	input: 'logs/{sample}_CELL_barcode.txt'
	output:
		pdf = 'plots/{sample}_CELL_dropped.pdf',
		png = 'plots/png/{sample}_CELL_dropped.png'
	params:
		min_quality = config['FILTER']['Cell_barcode']['min_quality'],
		num_below_quality = config['FILTER']['Cell_barcode']['num_below_quality']

	script:
		'../scripts/plot_cell_drop.R'

rule plot_BC_drop:
	input:
		BC_tagged = expand('logs/{sample}_CELL_barcode.txt', sample=samples.index),
		UMI_tagged = expand('logs/{sample}_UMI_barcode.txt', sample=samples.index),
		reads_left = expand('logs/{sample}_reads_left.txt', sample=samples.index)
	output:
		pdf = 'plots/BC_drop.pdf',
		png = 'plots/png/BC_drop.png'
	params:
		BC_length = config['FILTER']['Cell_barcode']['end'] - config['FILTER']['Cell_barcode']['start']+1,
		UMI_length = config['FILTER']['UMI']['end'] - config['FILTER']['UMI']['start']+1,
		min_num_below_BC = config['FILTER']['Cell_barcode']['num_below_quality'],
		min_num_below_UMI = config['FILTER']['UMI']['num_below_quality'],
		min_BC_quality = config['FILTER']['Cell_barcode']['min_quality'],
		min_UMI_quality = config['FILTER']['UMI']['min_quality'],
		sample_names = lambda wildcards: samples.index
	script:
		'../scripts/plot_BC_drop.R'

rule multiqc_trimmomatic:
	input: expand('logs/{sample}_trimlog.txt', sample=samples.index)
	output: 'reports/filter.html'
	shell:
		"""multiqc logs/ -m trimmomatic -o reports/ -n filter.html -f"""