"""Align the data with STAR."""


#Which rules will be run on the host computer and not sent to nodes
localrules: multiqc_star, plot_yield, plot_knee_plot, extend_barcode


rule STAR_align:
	input:
		fq1="data/{sample}/trimmmed_repaired_R2.fastq.gz",
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
		"0.27.1/bio/star/align"

# rule alevin:
# 	input:
# 		index='{salmon_index}',
# 		R1="data/{sample}/trimmmed_repaired_R1.fastq.gz",
# 		R2="data/{sample}/trimmmed_repaired_R2.fastq.gz",
# 	conda: '../envs/salmon.yaml'
# 	params:
# 		cell_barcode_length=(config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1),
# 		umi_barcode_length=(config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'] + 1)
# 	output:
# 		out_folder='data/{sample}/salmon/',
# 		counts='data/{sample}/salmon/mapping.tsv'
# 	shell:
# 		"""salmon alevin\
# 		-l ISR\
# 		-1 {input.R1}\
# 		-2 {input.R2}\
# 		-i {inout.index}\
# 		-p 10\
# 		-o {output.out_folder}\
# 		--tgMap {output.counts}\
# 		--barcodeLength {params.cell_barcode_length}\
# 		--umiLength {params.umi_barcode_length}\
# 		--end 5"""


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
		R1_ref = "data/{sample}/trimmmed_repaired_R1.fastq.gz"
	output:
		temp('data/{sample}/Aligned.merged.bam')
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
		data=temp('data/{sample}/Aligned.repaired.bam'),
		refFlat='{}.refFlat'.format(annotation_prefix)
	params:
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		temp('data/{sample}/gene_exon_tagged.bam')
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && TagReadWithGeneFunction -m {params.memory}\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}
		"""

rule DetectBeadSubstitutionErrors:
    input:
        'data/{sample}/gene_exon_tagged.bam'
    output:
        data=temp('data/{sample}/gene_exon_tagged_bead_sub.bam'),
        report='logs/{sample}_beadSubstitutionReport.txt',
        stats='logs/{sample}_beadSubstitutionStats.txt',
        summary='logs/{sample}_beadSubstitutionSummary.txt'
    params:
        SmartAdapter=config['FILTER']['5-prime-smart-adapter'],
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DetectBeadSynthesisErrors -m {params.memory}\
        I={input}\
        O={output.data}\
        REPORT={output.report}\
        OUTPUT_STATS={output.stats}\
        SUMMARY={output.summary}\
        PRIMER_SEQUENCE={params.SmartAdapter}
        """

rule bead_errors_metrics:
    input:
        'data/{sample}/gene_exon_tagged_bead_sub.bam'
    output:
        'data/{sample}/final.bam'
    params:
        out_stats='logs/{sample}_synthesis_stats.txt',
        summary='logs/{sample}_synthesis_stats_summary.txt',
        barcodes=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']) * 2,
        memory =config['LOCAL']['memory'],
        SmartAdapter=config['FILTER']['5-prime-smart-adapter'],
        temp_directory=config['LOCAL']['temp-directory']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DetectBeadSynthesisErrors -m {params.memory}\
        INPUT={input}\
        OUTPUT={output}\
        OUTPUT_STATS={params.out_stats}\
        SUMMARY={params.summary}\
        NUM_BARCODES={params.barcodes}\
        PRIMER_SEQUENCE={params.SmartAdapter}
        """


rule bam_hist:
	input:
		'data/{sample}/final.bam'
	params:
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'logs/dropseq_tools/{sample}_hist_out_cell.txt'
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && BamTagHistogram -m {params.memory}\
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
		data='logs/dropseq_tools/{sample}_hist_out_cell.txt',
		barcodes='data/{sample}/barcodes.csv'
	params: 
		cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells'])
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/knee_plots/{sample}_knee_plot.pdf'
	script:
		'../scripts/plot_knee_plot.R'
