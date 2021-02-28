"""Align the data with STAR."""


#Which rules will be run on the host computer and not sent to nodes
localrules:
    multiqc_star,
    plot_yield,
    pigz_unmapped,
    mv_outs_mtx,
    compress_mtx_out,
    plot_rna_metrics,
    SingleCellRnaSeqMetricsCollector,
    extract_bc_umi_data

ruleorder: STAR_solo_align_wl > STAR_solo_align

rule STAR_solo_align_wl:
    input:
        fq2='{results_dir}/samples/{sample}/trimmed_repaired_R2.fastq.gz',
        fq1='{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz',
        whitelist=config['FILTER']['barcode-whitelist'],
        index=lambda wildcards: '{}/{}_{}_{}/STAR_INDEXES/'.format(
            config['META']['reference-directory'],
            species,
            build,
            release),

    output:
        bam='{results_dir}/samples/{sample}/Aligned.sortedByCoord.out.bam',
        unmapped_R1='{results_dir}/samples/{sample}/Unmapped.out.mate1',
        unmapped_R2='{results_dir}/samples/{sample}/Unmapped.out.mate2',
        logs='{results_dir}/samples/{sample}/Log.final.out',
        solo_barcodes_raw='{results_dir}/samples/{sample}/Solo.out/Gene/raw/barcodes.tsv',
        solo_features_raw='{results_dir}/samples/{sample}/Solo.out/Gene/raw/features.tsv',
        solo_matrix_raw='{results_dir}/samples/{sample}/Solo.out/Gene/raw/matrix.mtx',
        solo_barcodes_filtered='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/barcodes.tsv',
        solo_features_filtered='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/features.tsv',
        solo_matrix_filtered='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/matrix.mtx',
    log:
        solo_barcode_stats = '{results_dir}/samples/{sample}/Solo.out/Barcodes.stats',
        solo_summary = '{results_dir}/samples/{sample}/Solo.out/Gene/Summary.csv',
        stdout = '{results_dir}/log/STARsolo/{sample}_stdout_stderr.log'
    params:
        extra="""--outSAMtype BAM SortedByCoordinate\
                --soloBarcodeReadLength 0\
                --outReadsUnmapped Fastx\
                --outFilterMismatchNmax {}\
                --outFilterMismatchNoverLmax {}\
                --outFilterMismatchNoverReadLmax {}\
                --outFilterMatchNmin {}\
                --outFilterScoreMinOverLread {}\
                --outFilterMatchNminOverLread {}\
                --outSAMattributes CR CY UR UY NH HI AS nM jM CB UB GX GN XS ch""".format(
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
            release),
        extra_solo=lambda wildcards:  """--soloFeatures Gene GeneFull\
            --soloType CB_UMI_Simple\
            --soloCBstart {}\
            --soloCBlen {}\
            --soloUMIstart {}\
            --soloUMIlen {}\
            --soloAdapterSequence {}\
            --soloCBmatchWLtype 1MM multi pseudocounts\
            --soloUMIdedup 1MM_Directional\
            --soloUMIfiltering MultiGeneUMI\
            --soloCellFilter TopCells {}""".format(
                config['FILTER']['cell-barcode']['start'],
                config['FILTER']['cell-barcode']['end'],
                config['FILTER']['UMI-barcode']['start'],
                1 + config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'],
                config['FILTER']['5-prime-smart-adapter'],
                samples.loc[wildcards.sample,'expected_cells']
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
            --soloCBwhitelist {input.whitelist}\
            2> {log.stdout}
      """

rule STAR_solo_align:
    input:
        fq2='{results_dir}/samples/{sample}/trimmed_repaired_R2.fastq.gz',
        fq1='{results_dir}/samples/{sample}/trimmed_repaired_R1.fastq.gz',
        index=lambda wildcards: '{}/{}_{}_{}/STAR_INDEXES/'.format(
            config['META']['reference-directory'],
            species,
            build,
            release),

    output:
        bam='{results_dir}/samples/{sample}/Aligned.sortedByCoord.out.bam',
        unmapped_R1='{results_dir}/samples/{sample}/Unmapped.out.mate1',
        unmapped_R2='{results_dir}/samples/{sample}/Unmapped.out.mate2',
        logs='{results_dir}/samples/{sample}/Log.final.out',
        solo_barcodes_raw='{results_dir}/samples/{sample}/Solo.out/Gene/raw/barcodes.tsv',
        solo_features_raw='{results_dir}/samples/{sample}/Solo.out/Gene/raw/features.tsv',
        solo_matrix_raw='{results_dir}/samples/{sample}/Solo.out/Gene/raw/matrix.mtx',
        solo_barcodes_filtered='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/barcodes.tsv',
        solo_features_filtered='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/features.tsv',
        solo_matrix_filtered='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/matrix.mtx',
    log:
        solo_barcode_stats = '{results_dir}/samples/{sample}/Solo.out/Barcodes.stats',
        solo_summary = '{results_dir}/samples/{sample}/Solo.out/Gene/Summary.csv',
        stdout = '{results_dir}/log/STARsolo/{sample}_stdout_stderr.log'
    params:
        extra="""--outSAMtype BAM SortedByCoordinate\
                --soloBarcodeReadLength 0\
                --outReadsUnmapped Fastx\
                --outFilterMismatchNmax {}\
                --outFilterMismatchNoverLmax {}\
                --outFilterMismatchNoverReadLmax {}\
                --outFilterMatchNmin {}\
                --outFilterScoreMinOverLread {}\
                --outFilterMatchNminOverLread {}\
                --outSAMattributes CR CY UR UY NH HI AS nM jM CB UB GX GN XS ch""".format(
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
            release),
        extra_solo=lambda wildcards:  """--soloFeatures Gene GeneFull\
            --soloType CB_UMI_Simple\
            --soloCBstart {}\
            --soloCBlen {}\
            --soloUMIstart {}\
            --soloUMIlen {}\
            --soloAdapterSequence {}\
            --soloCBmatchWLtype 1MM multi pseudocounts\
            --soloUMIdedup 1MM_Directional\
            --soloUMIfiltering MultiGeneUMI\
            --soloCellFilter TopCells {}""".format(
                config['FILTER']['cell-barcode']['start'],
                config['FILTER']['cell-barcode']['end'],
                config['FILTER']['UMI-barcode']['start'],
                1 + config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'],
                config['FILTER']['5-prime-smart-adapter'],
                samples.loc[wildcards.sample,'expected_cells']
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
            --soloCBwhitelist None\
            2> {log.stdout}
      """
    
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
        R1='{results_dir}/samples/{sample}/Unmapped.out.mate1',
        R2='{results_dir}/samples/{sample}/Unmapped.out.mate2'
    output:
        R1='{results_dir}/samples/{sample}/Unmapped.out.mate1.gz',
        R2='{results_dir}/samples/{sample}/Unmapped.out.mate2.gz'
    threads: 4
    conda: '../envs/pigz.yaml'
    shell:
        """pigz -p 4 {input.R1} {input.R2}"""

rule compress_mtx_out:
    input: 
        solo_barcodes='{results_dir}/samples/{sample}/Solo.out/Gene/raw/barcodes.tsv',
        solo_features='{results_dir}/samples/{sample}/Solo.out/Gene/raw/features.tsv',
        solo_matrix='{results_dir}/samples/{sample}/Solo.out/Gene/raw/matrix.mtx',
    output:
        solo_barcodes='{results_dir}/samples/{sample}/Solo.out/Gene/raw/barcodes.tsv.gz',
        solo_features='{results_dir}/samples/{sample}/Solo.out/Gene/raw/features.tsv.gz',
        solo_matrix='{results_dir}/samples/{sample}/Solo.out/Gene/raw/matrix.mtx.gz',
    conda: '../envs/pigz.yaml'
    threads: 3
    shell:
        """pigz -p {threads} {input.solo_barcodes} {input.solo_features} {input.solo_matrix}"""

rule mv_outs_mtx:
    input:
        solo_barcodes='{results_dir}/samples/{sample}/Solo.out/Gene/raw/barcodes.tsv.gz',
        solo_features='{results_dir}/samples/{sample}/Solo.out/Gene/raw/features.tsv.gz',
        solo_matrix='{results_dir}/samples/{sample}/Solo.out/Gene/raw/matrix.mtx.gz',
    output:
        '{results_dir}/samples/{sample}/umi/barcodes.tsv.gz',
        '{results_dir}/samples/{sample}/umi/features.tsv.gz',
        '{results_dir}/samples/{sample}/umi/matrix.mtx.gz',
    params:
        destination = "{results_dir}/samples/{sample}/umi"
    shell:
        """mv {input.solo_barcodes} {input.solo_features} {input.solo_matrix} {params.destination}"""


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


# rule plot_knee_plot:
#     input:
#         data='{results_dir}/logs/dropseq_tools/{sample}_hist_out_cell.txt',
#         barcodes='{results_dir}/samples/{sample}/filtered_barcodes.csv'
#     params:
#         cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells'])
#     conda: '../envs/r.yaml'
#     output:
#         pdf='{results_dir}/plots/knee_plots/{sample}_knee_plot.pdf'
#     script:
#         '../scripts/plot_knee_plot.R'

rule SingleCellRnaSeqMetricsCollector:
    input:
        data='{results_dir}/samples/{sample}/Aligned.sortedByCoord.out.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/barcodes.tsv',
        refFlat=expand("{ref_path}/{species}_{build}_{release}/curated_annotation.gtf",
            ref_path=ref_path,
            release=release,
            species=species,
            build=build),
        rRNA_intervals=expand("{ref_path}/{species}_{build}_{release}/annotation.rRNA.intervals",
            ref_path=ref_path,
            release=release,
            build=build,
            species=species)
    params:
        cells=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']),        
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        '{results_dir}/logs/dropseq_tools/{sample}_rna_metrics.txt'
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && SingleCellRnaSeqMetricsCollector -m {params.memory}\
        INPUT={input.data}\
        OUTPUT={output}\
        ANNOTATIONS_FILE={input.refFlat}\
        CELL_BC_FILE={input.barcode_whitelist}\
        RIBOSOMAL_INTERVALS={input.rRNA_intervals}
        """
rule plot_rna_metrics:
    input:
        rna_metrics='{results_dir}/logs/dropseq_tools/{sample}_rna_metrics.txt',
        barcode='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/barcodes.tsv'
    conda: '../envs/r.yaml'
    output:
        pdf='{results_dir}/plots/rna_metrics/{sample}_rna_metrics.pdf'
    script:
        '../scripts/plot_rna_metrics.R'

rule extract_bc_umi_data_all:
    input:
        expand('{results_dir}/samples/{sample}/bc_umi_long.csv', sample=samples.index, results_dir=results_dir)


rule extract_bc_umi_data:
    input: 
        data='{results_dir}/samples/{sample}/final.bam',
        barcode_whitelist='{results_dir}/samples/{sample}/Solo.out/Gene/filtered/barcodes.tsv'
    output:
        data='{results_dir}/samples/{sample}/bc_umi_long.csv'
    params:
        count_per_umi=config['EXTRACTION']['minimum-counts-per-UMI'],
        umiBarcodeEditDistance=config['EXTRACTION']['UMI-edit-distance'],
        temp_directory=config['LOCAL']['temp-directory'],
        memory=config['LOCAL']['memory'],
        locus_list=','.join(config['EXTRACTION']['LOCUS']),
        strand_strategy=config['EXTRACTION']['strand-strategy']
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && GatherMolecularBarcodeDistributionByGene -m {params.memory}\
        I={input.data}\
        O={output.data}\
        EDIT_DISTANCE={params.umiBarcodeEditDistance}\
        STRAND_STRATEGY={params.strand_strategy}\
        CELL_BARCODE_TAG=CR\
        MOLECULAR_BARCODE_TAG=UR\
        LOCUS_FUNCTION_LIST=null\
        LOCUS_FUNCTION_LIST={{{params.locus_list}}}\
        MIN_BC_READ_THRESHOLD={params.count_per_umi}\
        CELL_BC_FILE={input.barcode_whitelist}"""