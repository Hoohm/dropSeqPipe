"""Align the data with STAR."""

ruleorder: plot_knee_plot_whitelist > plot_knee_plot


#Which rules will be run on the host computer and not sent to nodes
localrules: multiqc_star, plot_yield, plot_knee_plot, plot_knee_plot_whitelist


rule STAR_align:
	input:
		fq1="data/{sample}/{sample}_trimmmed_repaired_R2.fastq.gz",
		index=lambda wildcards: star_index_prefix + '_' + str(samples.loc[wildcards.sample,'read_length']) + '/SA'
	output:
		temp('data/{sample}/Aligned.out.bam')
	log:
		'data/{sample}/Log.final.out'
	params:
		extra="""--outReadsUnmapped Fastx\
				--outSAMattributes NH HI NM AS MD\
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
		mapped='data/{sample}/Aligned.out.bam',
		R1_ref = "data/{sample}/{sample}_trimmmed_repaired_R1.fastq.gz"
	output:
		temp('data/{sample}.Aligned.merged.bam')
	params:
		BC_start=config['FILTER']['cell-barcode']['start']-1,
		BC_end=config['FILTER']['cell-barcode']['end'],
		UMI_start=config['FILTER']['UMI-barcode']['start']-1,
		UMI_end=config['FILTER']['UMI-barcode']['end'],
		discard_secondary_alignements=True
	conda: '../envs/merge_bam.yaml'
	script:
		'../scripts/merge_bam.py'

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
		R1_filtered=expand('logs/cutadapt/{sample}_R1.qc.txt', sample=samples.index),
		R2_filtered=expand('logs/cutadapt/{sample}_R2.qc.txt', sample=samples.index),
		repaired=expand('logs/bbmap/{sample}_repair.txt', sample=samples.index),
		STAR_output=expand('data/{sample}/Log.final.out', sample=samples.index),
	params:
		BC_length=config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start']+1,
		UMI_length=config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start']+1,
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

rule violine_plots:
	input:
		UMIs='summary/umi_expression_matrix.tsv',
		counts='summary/counts_expression_matrix.tsv',
		design='samples.csv'
#	params:
#		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
#		edit_distance=config['EXTRACTION']['UMI-edit-distance']
	conda: '../envs/plots_ext.yaml'
	output:
		pdf_violine='plots/violinplots_comparison_UMI.pdf',
#		html_umivscounts='plots/UMI_vs_counts.html',
		pdf_umivscounts='plots/UMI_vs_counts.pdf',
#		html_umi_vs_gene='plots/UMI_vs_gene.html',
		pdf_umi_vs_gene='plots/UMI_vs_gene.pdf',
#		html_count_vs_gene='plots/Count_vs_gene.html',
		pdf_count_vs_gene='plots/Count_vs_gene.pdf',
		R_objects='summary/R_Seurat_objects.rdata'
	script:
		'../scripts/plot_violine.R'
