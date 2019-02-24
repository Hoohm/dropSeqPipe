import math
import platform
"""Generate all the meta data files"""
# To add missing fields for an annotation of ERCC: awk -F'[\t|;]' '{printf $0" "; gsub(/id/,"name"); print $9";"$10"; exon_version \"1\";"}'
#Which rules will be run on the host computer and not sent to nodes
localrules:
     create_dict,
     reduce_gtf,
     create_refFlat,
     create_intervals,
     curate_annotation


rule curate_annotation:
    input:
        biotypes=config['META']['gtf_biotypes'],
        annotation="{ref_path}/{species}_{build}_{release}/annotation.gtf"
    output:
        temp("{ref_path}/{species}_{build}_{release}/curated_annotation.gtf")
    params:
        patterns='|'.join(config['biotypes'])
    shell:
        """cat {input.annotation} | grep -E "{params.patterns}" > {output}"""


rule create_dict:
    input:
        "{ref_path}/{species}_{build}_{release}/genome.fa"
    output:
        "{ref_path}/{species}_{build}_{release}/genome.dict"
    threads:1
    params:
        picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
        temp_directory=config['LOCAL']['temp-directory']
    conda: '../envs/picard.yaml'
    shell:
        """java -jar -Djava.io.tmpdir={params.temp_directory} {params.picard} CreateSequenceDictionary\
        REFERENCE={input}\
        OUTPUT={output}
        """

rule reduce_gtf:
    input:
        reference_dict="{ref_path}/{species}_{build}_{release}/genome.dict",
        annotation="{ref_path}/{species}_{build}_{release}/curated_annotation.gtf"
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        "{ref_path}/{species}_{build}_{release}/curated_reduced_annotation.gtf"
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ReduceGtf -m {params.memory}\
        GTF={input.annotation}\
        OUTPUT={output}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        IGNORE_FUNC_TYPE='null'\
        ENHANCE_GTF='false'"""

rule create_refFlat:
    input:
        reference_dict="{ref_path}/{species}_{build}_{release}/genome.dict",
        annotation="{ref_path}/{species}_{build}_{release}/curated_annotation.gtf"
    params:
        memory=config['LOCAL']['memory'],
        temp_directory=config['LOCAL']['temp-directory']
    output:
        "{ref_path}/{species}_{build}_{release}/curated_annotation.refFlat"
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ConvertToRefFlat -m {params.memory}\
        ANNOTATIONS_FILE={input.annotation}\
        OUTPUT={output}\
        SEQUENCE_DICTIONARY={input.reference_dict}
        """

rule create_intervals:
    input:
        annotation_reduced="{ref_path}/{species}_{build}_{release}/curated_reduced_annotation.gtf",
        reference_dict="{ref_path}/{species}_{build}_{release}/genome.dict"
    params:
        memory=config['LOCAL']['memory'],
        reference_directory=config['META']['reference-directory'],
        temp_directory=config['LOCAL']['temp-directory'],
        prefix="{species}_{build}_{release}/annotation"
    output:
        intervals="{ref_path}/{species}_{build}_{release}/annotation.rRNA.intervals"
    conda: '../envs/dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && CreateIntervalsFiles -m {params.memory}\
        REDUCED_GTF={input.annotation_reduced}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        O={params.reference_directory}\
        PREFIX={params.prefix}
        """

rule get_genomeChrBinNbits:
    input:
        reference_file="{ref_path}/{species}_genome.fa"
    params:
        samples_file='samples.csv',
        reference_directory=config['META']['reference-directory']
    output:
        '{params.reference_directory}/index_params.txt'
    run:
        """
        from math import log2
        from platform import system
        if (system() == 'Darwin'):
            genomeLength = shell("wc -c {} | cut -d' ' -f2".format(snakemake.reference_file), iterable=True)
        else:
            genomeLength = shell("wc -c {} | cut -d' ' -f1".format(snakemake.reference_file), iterable=True)
        genomeLength = int(next(genomeLength))
        referenceNumber = shell('grep "^>" {} | wc -l'.format(snakemake.reference_file), iterable=True)
        referenceNumber = int(next(referenceNumber))
        value = min([18,int(log2(genomeLength/referenceNumber))])
        """

def get_sjdbOverhang(wildcards):
    return(int(wildcards.read_length)-1)


rule prep_star_index:
    input:
        reference_file="{ref_path}/{species}_genome.fa",
        config_file='config.yaml'
    output:
        '{reference_directory}/star_ref_config.txt'
    conda:
        '../envs/pyyaml.yaml'
    script:
        '../scripts/prep_star.py'


    

rule create_star_index:
    input:
        reference_file="{ref_path}/{species}_{build}_{release}/genome.fa",
        annotation_file="{ref_path}/{species}_{build}_{release}/curated_annotation.gtf"  
    params:
        sjdbOverhang=lambda wildcards: get_sjdbOverhang(wildcards),
        genomeDir='{ref_path}/{species}_{build}_{release}/STAR_INDEX/SA_{read_length}',
        genomeChrBinNbits=config['MAPPING']['STAR']['genomeChrBinNbits']
    output:
        '{ref_path}/{species}_{build}_{release}/STAR_INDEX/SA_{read_length}/SA'
    threads: 24
    conda: '../envs/star.yaml'
    shell:
        """mkdir -p {params.genomeDir}; STAR\
        --runThreadN {threads}\
        --runMode genomeGenerate\
        --genomeDir {params.genomeDir}\
        --genomeFastaFiles {input.reference_file}\
        --sjdbGTFfile {input.annotation_file}\
        --sjdbOverhang {params.sjdbOverhang}\
        --genomeChrBinNbits {params.genomeChrBinNbits}\
        --genomeSAsparseD 2
        """