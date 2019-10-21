"""Align the data with STAR."""


#Which rules will be run on the host computer and not sent to nodes
localrules:
    multiqc_star,
    plot_yield,
    plot_knee_plot,
    pigz_unmapped


rule STAR_align:
    input:
        fq1='{results_dir}/samples/{sample}/trimmed_repaired_R2.fastq.gz',
        index=lambda wildcards: '{}/{}_{}_{}/STAR_INDEX/SA'.format(
            config['META']['reference-directory'],
            species,
            build,
            release) + '_' + str(samples.loc[wildcards.sample,'read_length']) + '/SA'
    output:
        temp('{results_dir}/samples/{sample}/Aligned.out.bam'),
        '{results_dir}/samples/{sample}/Unmapped.out.mate1'

    log:
        '{results_dir}/samples/{sample}/Log.final.out'
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
        index=lambda wildcards: '{}/{}_{}_{}/STAR_INDEX/SA'.format(
            config['META']['reference-directory'],
            species,
            build,
            release) + '_' + str(samples.loc[wildcards.sample,'read_length']) + '/'
    singularity:
        "shub://seb-mueller/singularity_dropSeqPipe:v04"
    threads: 24
    wrapper:
        "0.27.1/bio/star/align"
# rule alevin:
#   input:
#       index='{salmon_index}',
#       R1="samples/{sample}/trimmed_repaired_R1.fastq.gz",
#       R2="samples/{sample}/trimmed_repaired_R2.fastq.gz",
#   conda: '../envs/salmon.yaml'
#   params:
#       cell_barcode_length=(config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1),
#       umi_barcode_length=(config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'] + 1)
#   output:
#       out_folder='samples/{sample}/salmon/',
#       counts='samples/{sample}/salmon/mapping.tsv'
#   shell:
#       """salmon alevin\
#       -l ISR\
#       -1 {input.R1}\
#       -2 {input.R2}\
#       -i {inout.index}\
#       -p 10\
#       -o {output.out_folder}\
#       --tgMap {output.counts}\
#       --barcodeLength {params.cell_barcode_length}\
#       --umiLength {params.umi_barcode_length}\
#       --end 5"""


rule multiqc_star:
    input:
        expand('{results_dir}/samples/{sample}/Log.final.out', sample=samples.index, results_dir=results_dir)
    output:
        html='{results_dir}/reports/star.html'
    params: '-m star'
    wrapper:
        '0.36.0/bio/multiqc'

rule pigz_unmapped:
    input:
        '{results_dir}/samples/{sample}/Unmapped.out.mate1'
    output:
        '{results_dir}/samples/{sample}/Unmapped.out.mate1.gz'
    threads: 4
    conda: '../envs/pigz.yaml'
    shell:
        """pigz -p 4 {input}"""

rule MergeBamAlignment:
    input:
        mapped='{results_dir}/samples/{sample}/Aligned.out.bam',
        R1_ref = '{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz'
    output:
        temp('{results_dir}/samples/{sample}/Aligned.merged.bam')
    params:
        BC_start=config['FILTER']['cell-barcode']['start']-1,
        BC_end=config['FILTER']['cell-barcode']['end'],
        UMI_start=config['FILTER']['UMI-barcode']['start']-1,
        UMI_end=config['FILTER']['UMI-barcode']['end'],
        discard_secondary_alignements=True
    conda: '../envs/merge_bam.yaml'
    script:
        '../scripts/merge_bam.py'

# Note: rule repair_barcodes (cell_barcodes.smk) creates Aligned.repaired.bam
# this is using barcode information (i.e. dependent on expected_cells in config.yaml)


rule TagReadWithGeneExon:
    input:
        data='{results_dir}/samples/{sample}/Aligned.repaired.bam',
        refFlat=expand("{ref_path}/{species}_{build}_{release}/curated_annotation.refFlat",
            ref_path=config['META']['reference-directory'],
            species=species,
            release=release,
            build=build)
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        temp('{results_dir}/samples/{sample}/gene_exon_tagged.bam')
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && TagReadWithGeneFunction -m {params.memory}\
        INPUT={input.data}\
        OUTPUT={output}\
        ANNOTATIONS_FILE={input.refFlat}
        """

rule DetectBeadSubstitutionErrors:
    input:
        '{results_dir}/samples/{sample}/gene_exon_tagged.bam'
    output:
        data=temp('{results_dir}/samples/{sample}/gene_exon_tagged_bead_sub.bam'),
        report='{results_dir}/logs/dropseq_tools/{sample}_beadSubstitutionReport.txt',
        summary='{results_dir}/logs/dropseq_tools/{sample}_beadSubstitutionSummary.txt'
    params:
        SmartAdapter=config['FILTER']['5-prime-smart-adapter'],
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    conda: '../envs/dropseq_tools.yaml'
    threads: 5
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DetectBeadSubstitutionErrors -m {params.memory}\
        I={input}\
        O={output.data}\
        OUTPUT_REPORT={output.report}\
        OUTPUT_SUMMARY={output.summary}\
        NUM_THREADS={threads}
        """

rule bead_errors_metrics:
    input:
        '{results_dir}/samples/{sample}/gene_exon_tagged_bead_sub.bam'
    output:
        '{results_dir}/samples/{sample}/final.bam'
    params:
        out_stats='{results_dir}/logs/dropseq_tools/{sample}_synthesis_stats.txt',
        summary='{results_dir}/logs/dropseq_tools/{sample}_synthesis_stats_summary.txt',
        barcodes=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']) * 2,
        memory =config['LOCAL']['memory'],
        SmartAdapter=config['FILTER']['5-prime-smart-adapter'],
        temp_directory=config['LOCAL']['temp-directory']
    conda: '../envs/dropseq_tools.yaml'
    threads: 5
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DetectBeadSynthesisErrors -m {params.memory}\
        INPUT={input}\
        OUTPUT={output}\
        OUTPUT_STATS={params.out_stats}\
        SUMMARY={params.summary}\
        NUM_BARCODES={params.barcodes}\
        PRIMER_SEQUENCE={params.SmartAdapter}\
        NUM_THREADS={threads}
        """


rule bam_hist:
    input:
        '{results_dir}/samples/{sample}/final.bam'
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        '{results_dir}/logs/dropseq_tools/{sample}_hist_out_cell.txt'
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && BamTagHistogram -m {params.memory}\
        TAG=XC\
        I={input}\
        READ_MQ=10\
        O={output}
        """


rule plot_yield:
    input:
        R1_filtered=expand('{results_dir}/logs/cutadapt/{sample}_R1.qc.txt', sample=samples.index, results_dir=results_dir),
        R2_filtered=expand('{results_dir}/logs/cutadapt/{sample}_R2.qc.txt', sample=samples.index, results_dir=results_dir),
        repaired=expand('{results_dir}/logs/bbmap/{sample}_repair.txt', sample=samples.index, results_dir=results_dir),
        STAR_output=expand('{results_dir}/samples/{sample}/Log.final.out', sample=samples.index, results_dir=results_dir),
    params:
        BC_length=config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start']+1,
        UMI_length=config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start']+1,
        sample_names=lambda wildcards: samples.index,
        batches=lambda wildcards: samples.loc[samples.index, 'batch']
    conda: '../envs/r.yaml'
    output:
        pdf='{results_dir}/plots/yield.pdf'
    script:
        '../scripts/plot_yield.R'


rule plot_knee_plot:
    input:
        data='{results_dir}/logs/dropseq_tools/{sample}_hist_out_cell.txt',
        barcodes='{results_dir}/samples/{sample}/barcodes.csv'
    params:
        cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells'])
    conda: '../envs/r.yaml'
    output:
        pdf='{results_dir}/plots/knee_plots/{sample}_knee_plot.pdf'
    script:
        '../scripts/plot_knee_plot.R'
