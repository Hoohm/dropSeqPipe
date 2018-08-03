"""Align the data with STAR."""

ruleorder: plot_knee_plot_whitelist > plot_knee_plot


#Which rules will be run on the host computer and not sent to nodes
localrules: multiqc_star, plot_yield, plot_knee_plot, plot_knee_plot_whitelist


rule STAR_align:
	input:
		fq1="data/{sample}_filtered.fastq.gz",
		index=lambda wildcards: star_index_prefix + '_' + str(samples.loc[wildcards.sample,'read_length']) + '/SA'
	output:
		temp('data/{sample}/Aligned.out.bam')
	log:
		'data/{sample}/Log.final.out'
	params:
		extra="""--outReadsUnmapped Fastx\
			 	--outFilterMismatchNmax {}\
			 	--outFilterMismatchNoverLmax {}\
			 	--outFilterMismatchNoverReadLmax {}\
			 	--outFilterMatchNmin {}\
			 	--outFilterScoreMinOverLread {}\
			 	--outFilterMatchNminOverLread {}""".format(
				config['MAPPING']['STAR']['outFilterMismatchNmax'],
				config['MAPPING']['STAR']['outFilterMismatchNoverLmax'],
				config['MAPPING']['STAR']['outFilterMismatchNoverReadLmax'],
				config['MAPPING']['STAR']['outFilterMatchNmin'],
				config['MAPPING']['STAR']['outFilterMatchNminOverLread'],
				config['MAPPING']['STAR']['outFilterScoreMinOverLread'],),
		index=lambda wildcards: star_index_prefix + '_' + str(samples.loc[wildcards.sample,'read_length']) + '/'	
	threads: 24
	wrapper:
		"0.22.0/bio/star/align"

rule sort_sam:
	input:
		'data/{sample}/Aligned.out.bam'
	params:
		temp_directory=config['LOCAL']['temp-directory'],
		picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
		memory=config['LOCAL']['memory']
	output:
		temp('data/{sample}.Aligned.sorted.bam')
	conda: '../envs/picard.yaml'
	shell:
		"""java -Xmx{params.memory} -jar -Djava.io.tmpdir={params.temp_directory} {params.picard} SortSam\
		INPUT={input}\
		OUTPUT={output}\
		SORT_ORDER=queryname\
		TMP_DIR={params.temp_directory}
		"""
rule multiqc_star:
	input:
		expand('data/{sample}/Log.final.out', sample=samples.index)
	output:
		html='reports/star.html'
	params: '-m star'
	wrapper:
		'0.21.0/bio/multiqc'


rule MergeBamAlignment:
	input:
		unmapped='data/{sample}_trimmed_unmapped.bam',
		mapped='data/{sample}.Aligned.sorted.bam',
		dict_file='{}.dict'.format(reference_prefix)
	output:
		temp('data/{sample}.Aligned.merged.bam')
	params:
		picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
		reference_file=reference_file,
		temp_directory=config['LOCAL']['temp-directory'],
		memory=config['LOCAL']['memory']
	conda: '../envs/picard.yaml'
	shell:
		"""java -Djava.io.tmpdir={params.temp_directory} -Xmx{params.memory} -jar {params.picard} MergeBamAlignment\
		REFERENCE_SEQUENCE={params.reference_file}\
		UNMAPPED_BAM={input.unmapped}\
		ALIGNED_BAM={input.mapped}\
		INCLUDE_SECONDARY_ALIGNMENTS=false\
		PAIRED_RUN=false\
		OUTPUT={output}
		"""
rule TagReadWithGeneExon:
	input:
		data='data/{sample}.Aligned.merged.bam',
		refFlat='{}.refFlat'.format(annotation_prefix)
	params:
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		temp('data/{sample}_gene_exon_tagged.bam')
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS="-Djava.io.tmpdir={params.temp_directory} $_JAVA_OPTIONS" && TagReadWithGeneExon -m {params.memory}\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		TAG=GE\
		CREATE_INDEX=true
		"""

rule bead_errors_metrics:
	input:
		'data/{sample}_gene_exon_tagged.bam'
	output:
		'data/{sample}_final.bam'
	params:
		out_stats='logs/{sample}_synthesis_stats.txt',
		summary='logs/{sample}_synthesis_stats_summary.txt',
		barcodes=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']) * 2,
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory =config['LOCAL']['memory'],
		SmartAdapter=config['FILTER']['5-prime-smart-adapter'],
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS="-Djava.io.tmpdir={params.temp_directory} $_JAVA_OPTIONS" && DetectBeadSynthesisErrors -m {params.memory}\
		INPUT={input}\
		OUTPUT={output}\
		OUTPUT_STATS={params.out_stats}\
		SUMMARY={params.summary}\
		NUM_BARCODES={params.barcodes}\
		PRIMER_SEQUENCE={params.SmartAdapter}
		"""

rule bam_hist:
	input:
		'data/{sample}_final.bam'
	params:
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'logs/{sample}_hist_out_cell.txt'
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS="-Djava.io.tmpdir={params.temp_directory} $_JAVA_OPTIONS" && BAMTagHistogram -m {params.memory}\
		TAG=XC\
		I={input}\
		READ_QUALITY=10\
		O={output}
		"""

rule plot_yield:
	input:
		BC_tagged=expand('logs/{sample}_CELL_barcode.txt', sample=samples.index),
		UMI_tagged=expand('logs/{sample}_UMI_barcode.txt', sample=samples.index),
		reads_left=expand('logs/{sample}_reads_left.txt', sample=samples.index),
		STAR_output=expand('data/{sample}/Log.final.out', sample=samples.index),
		trimmomatic_filtered=expand('logs/{sample}_reads_left_trim.txt', sample=samples.index)
	params:
		BC_length=config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start']+1,
		UMI_length=config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start']+1,
		min_num_below_BC=config['FILTER']['cell-barcode']['num-below-quality'],
		min_num_below_UMI=config['FILTER']['UMI-barcode']['num-below-quality'],
		min_BC_quality=config['FILTER']['cell-barcode']['min-quality'],
		min_UMI_quality=config['FILTER']['UMI-barcode']['min-quality'],
		sample_names=lambda wildcards: samples.index,
		batches=lambda wildcards: samples.loc[samples.index, 'batch']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/yield.pdf'
	script:
		'../scripts/plot_yield.R'

rule plot_knee_plot:
	input:
		'logs/{sample}_hist_out_cell.txt'
	params: 
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		edit_distance=config['EXTRACTION']['UMI-edit-distance']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/{sample}_knee_plot.pdf'
	script:
		'../scripts/plot_knee_plot.R'

rule plot_knee_plot_whitelist:
	input:
		data='logs/{sample}_hist_out_cell.txt',
		barcodes='barcodes.csv'
	params: 
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/{sample}_knee_plot.pdf'
	script:
		'../scripts/plot_knee_plot.R'