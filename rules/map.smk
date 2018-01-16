"""Align the data with STAR."""

ruleorder: plot_knee_plot_whitelist > plot_knee_plot

rule STAR_align:
	input: 
		data="data/{sample}_filtered.fastq.gz",
		index = expand('{star_index_prefix}_{read_length}/SA', star_index_prefix=star_index_prefix, read_length=read_lengths)
	output:
		sam = temp('logs/{sample}.Aligned.out.bam'),
		log_out = 'logs/{sample}.Log.final.out'
	params:
		prefix = '{sample}.',
		outFilterMismatchNmax = config['STAR_PARAMETERS']['outFilterMismatchNmax'],
		outFilterMismatchNoverLmax = config['STAR_PARAMETERS']['outFilterMismatchNoverLmax'],
		outFilterMismatchNoverReadLmax = config['STAR_PARAMETERS']['outFilterMismatchNoverReadLmax'],
		outFilterMatchNmin = config['STAR_PARAMETERS']['outFilterMatchNmin'],
		index = expand('{star_index_prefix}_{read_length}/', star_index_prefix=star_index_prefix, read_length=read_lengths)
	threads: 24
	shell:"""STAR\
			--genomeDir {params.index}\
			--readFilesCommand zcat\
			--runThreadN {threads}\
			--readFilesIn {input.data}\
			--outSAMtype BAM Unsorted\
			--outReadsUnmapped Fatsx\
			--outFileNamePrefix logs/{params.prefix}\
			--outFilterMismatchNmax {params.outFilterMismatchNmax}\
			--outFilterMismatchNoverLmax {params.outFilterMismatchNoverLmax}\
			--outFilterMismatchNoverReadLmax {params.outFilterMismatchNoverReadLmax}\
			--outFilterMatchNmin {params.outFilterMatchNmin}
			"""

rule sort_sam:
	input: 'logs/{sample}.Aligned.out.bam'
	params:
		TMPDIR = config['LOCAL']['TMPDIR'],
		picard = "$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
		memory = config['LOCAL']['MEMORY']
	output: temp('data/{sample}.Aligned.sorted.bam')
	shell:
		"""java -Xmx{params.memory} -jar -Djava.io.tmpdir={params.TMPDIR}	{params.picard} SortSam\
		INPUT={input}\
		OUTPUT={output}\
		SORT_ORDER=queryname\
		TMP_DIR={params.TMPDIR}
		"""

rule multiqc_star:
	input: expand('logs/{sample}.Log.final.out', sample = samples.index)
	output: 'reports/star.html'
	shell:
		"""multiqc logs/ -m star -o reports/ -n star.html -f"""


rule MergeBamAlignment:
	input:
		unmapped = 'data/{sample}_trimmed_unmapped.bam',
		mapped = 'data/{sample}.Aligned.sorted.bam',
		dict_file = expand('{reference_prefix}.dict', reference_prefix=reference_prefix)
	output: temp('data/{sample}.Aligned.merged.bam')
	params:
		picard = "$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
		reference = reference_file,
		TMPDIR = config['LOCAL']['TMPDIR'],
		memory = config['LOCAL']['MEMORY']
	threads: 1
	shell:
		"""java -Djava.io.tmpdir={params.TMPDIR} -Xmx{params.memory} -jar {params.picard} MergeBamAlignment\
		REFERENCE_SEQUENCE={params.reference}\
		UNMAPPED_BAM={input.unmapped}\
		ALIGNED_BAM={input.mapped}\
		INCLUDE_SECONDARY_ALIGNMENTS=false\
		PAIRED_RUN=false\
		OUTPUT={output}
		"""
rule TagReadWithGeneExon:
	input:
		data = 'data/{sample}.Aligned.merged.bam',
		refFlat = expand('{annotation_prefix}.refFlat', annotation_prefix=annotation_prefix)
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR']
	output:
		temp('data/{sample}_gene_exon_tagged.bam')
	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p TagReadWithGeneExon\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		TAG=GE\
		CREATE_INDEX=true
		"""

rule bead_errors_metrics:
	input: 'data/{sample}_gene_exon_tagged.bam'
	output:
		'data/{sample}_final.bam'
	params:
		out_stats = 'logs/{sample}_synthesis_stats.txt',
		summary = 'logs/{sample}_synthesis_stats_summary.txt',
		barcodes = lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']) * 2,
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory =config['LOCAL']['MEMORY'],
		SmartAdapter = config['FILTER']['5PrimeSmartAdapter'],
		TMPDIR = config['LOCAL']['TMPDIR']

	shell:
		"""{params.DROPSEQ_wrapper} -t {params.TMPDIR} -m {params.memory} -p DetectBeadSynthesisErrors\
		INPUT={input}\
		OUTPUT={output}\
		OUTPUT_STATS={params.out_stats}\
		SUMMARY={params.summary}\
		NUM_BARCODES={params.barcodes}\
		PRIMER_SEQUENCE={params.SmartAdapter}
		"""

rule bam_hist:
	input: 'data/{sample}_final.bam'
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR=config['LOCAL']['TMPDIR']
	output: 'logs/{sample}_hist_out_cell.txt'
	shell:
		"""{params.DROPSEQ_wrapper} -p BAMTagHistogram -m {params.memory} -t {params.TMPDIR}\
		TAG=XC\
		I={input}\
		O={output}
		"""

rule plot_knee_plot:
	input: 'logs/{sample}_hist_out_cell.txt'
	output:
		plot = 'plots/{sample}_knee_plot.pdf'
	params: 
		cells = lambda wildcards: samples.loc[wildcards.sample,'expected_cells']
	script:
		'../scripts/plot_knee_plot.R'

rule plot_knee_plot_whitelist:
	input:
		data = 'logs/{sample}_hist_out_cell.txt',
		barcodes = 'barcodes.csv'
	output:
		plot = 'plots/{sample}_knee_plot.pdf'
	params: 
		cells = lambda wildcards: samples.loc[wildcards.sample,'expected_cells']
	script:
		'../scripts/plot_knee_plot.R'