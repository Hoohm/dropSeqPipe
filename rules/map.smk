"""Align the data with STAR."""


#Which rules will be run on the host computer and not sent to nodes
localrules:
    multiqc_star,
    plot_yield,
    plot_knee_plot,
    pigz_unmapped


rule STAR_solo_align:
    input:
        fq2='{results_dir}/samples/{sample}/trimmed_repaired_R2.fastq.gz',
        fq1='{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz',
        index=lambda wildcards: '{}/{}_{}_{}/STAR_INDEXES/'.format(
            config['META']['reference-directory'],
            species,
            build,
            release) + str(samples.loc[wildcards.sample,'read_length']),
            whitelist='{results_dir}/samples/{sample}/barcodes.csv'

    output:
        bam=temp('{results_dir}/samples/{sample}/Aligned.sortedByCoord.out.bam'),
        unmapped='{results_dir}/samples/{sample}/Unmapped.out.mate1',
        logs='{results_dir}/samples/{sample}/Log.final.out'
    params:
        extra="""--outSAMtype BAM SortedByCoordinate\
                --outReadsUnmapped Fastx\
                --outFilterMismatchNmax {}\
                --outFilterMismatchNoverLmax {}\
                --outFilterMismatchNoverReadLmax {}\
                --outFilterMatchNmin {}\
                --outFilterScoreMinOverLread {}\
                --outFilterMatchNminOverLread {}\
                --outSAMattributes CR CY UR UY NH HI AS nM jM""".format(
                config['MAPPING']['STAR']['outFilterMismatchNmax'],
                config['MAPPING']['STAR']['outFilterMismatchNoverLmax'],
                config['MAPPING']['STAR']['outFilterMismatchNoverReadLmax'],
                config['MAPPING']['STAR']['outFilterMatchNmin'],
                config['MAPPING']['STAR']['outFilterMatchNminOverLread'],
                config['MAPPING']['STAR']['outFilterScoreMinOverLread'],),
        index=lambda wildcards: '{}/{}_{}_{}/STAR_INDEXES/'.format(
            config['META']['reference-directory'],
            species,
            build,
            release) + str(samples.loc[wildcards.sample,'read_length']),
        extra_solo="""--soloFeatures Gene\
            --soloType Droplet\
            --soloCBstart {}\
            --soloCBlen {}\
            --soloUMIstart {}\
            --soloUMIlen {}\
            --soloAdapterSequence {}\
            --soloCBmatchWLtype 1MM\
            --soloUMIdedup 1MM_Directional\
            --soloCellFilter None""".format(
                config['FILTER']['cell-barcode']['start'],
                config['FILTER']['cell-barcode']['end'],
                config['FILTER']['UMI-barcode']['start'],
                1 + config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'],
                config['FILTER']['5-prime-smart-adapter']
            ),
        out_prefix=lambda wildcards: '{}/samples/{}/'.format(wildcards.results_dir, wildcards.sample)
    singularity:
        "shub://seb-mueller/singularity_dropSeqPipe:v04"
    threads: 24
    conda: '../envs/star.yaml'
    shell:
        """STAR {params.extra} {params.extra_solo}\
            --runThreadN {threads}\
            --genomeDir {input.index}\
            --readFilesCommand zcat\
            --readFilesIn {input.fq2} {input.fq1}\
            --outFileNamePrefix {params.out_prefix}\
            --soloCBwhitelist {input.whitelist}"""

    
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


rule TagReadWithGeneExon:
    input:
        data='{results_dir}/samples/{sample}/Aligned.sortedByCoord.out.bam',
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
        CELL_BARCODE_TAG=CR\
        MOLECULAR_BARCODE_TAG=UR\
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
        CELL_BARCODE_TAG=CR\
        MOLECULAR_BARCODE_TAG=UR\
        NUM_THREADS={threads}
        """



rule plot_yield:
    input:
        R1_filtered=expand('{results_dir}/logs/cutadapt/{sample}_R1.qc.txt', sample=samples.index, results_dir=results_dir),
        R2_filtered=expand('{results_dir}/logs/cutadapt/{sample}_R2.qc.txt', sample=samples.index, results_dir=results_dir),
        repaired=expand('{results_dir}/logs/repair/{sample}.csv', sample=samples.index, results_dir=results_dir),
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
        barcodes='{results_dir}/samples/{sample}/filtered_barcodes.csv'
    params:
        cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells'])
    conda: '../envs/r.yaml'
    output:
        pdf='{results_dir}/plots/knee_plots/{sample}_knee_plot.pdf'
    script:
        '../scripts/plot_knee_plot.R'
