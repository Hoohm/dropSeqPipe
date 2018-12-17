from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

localrules:
    download_annotation,
    download_genome,
    rename_genome,
    merge_genomes,
    merge_annotations

def get_annotation(wildcards):
    return FTP.remote("ftp.ensembl.org/pub/release-{0}/gtf/{1}/{2}.GRC{3}{4}.{0}.gtf.gz".format(
        wildcards.release,
        wildcards.species,
        wildcards.species.capitalize(),
        wildcards.species.lower()[0],
        wildcards.build),
        static=True,
        keep_local=True)

def get_genome(wildcards):
    return FTP.remote("ftp.ensembl.org/pub/release-{0}/fasta/{1}/dna/{2}.GRC{3}{4}.dna.primary_assembly.fa.gz".format(
        wildcards.release,
        wildcards.species,
        wildcards.species.capitalize(),
        wildcards.species.lower()[0],
        wildcards.build),
        static=True,
        keep_local=True)

rule download_annotation:
    input:
        get_annotation
    output:
        "{ref_path}/{species}_{build}_{release}/annotation.gtf"
    shell:
        "gunzip -c -d {input} > {output}"

rule download_genome:
    input:
        get_genome
    output:
        "{ref_path}/{species}_{build}_{release}/genome.fa"
    shell:
        "gunzip -d -c {input} > {output}"


rule rename_genome:
    input:
        "{ref_path}/{species}_{build}_{release}/genome.fa"
    output:
       temp("{ref_path}/{species}_{build}_{release}/renamed_genome.fa")
    params:
        species= lambda wildcards: wildcards.species
    shell:
        """sed -e 's/>/>{params.species}/g' {input} > {output}"""


rule merge_genomes:
    input:
        genome1=expand("{ref_path}/{species}_{build}_{release}/renamed_genome.fa",
            species=species_list[0],
            build=build_list[0],
            release=release_list[0],
            ref_path=config['META']['reference-directory']),
        genome2=expand("{ref_path}/{species}_{build}_{release}/renamed_genome.fa",
            species=species_list[1],
            build=build_list[1],
            release=release_list[1],
            ref_path=config['META']['reference-directory'])
    output:
        "{}/{}_{}_{}/genome.fa".format(
            config['META']['reference-directory'],
            species,
            build,
            release)
    shell:
        """cat {input.genome1} {input.genome2} > {output}"""

rule merge_annotations:
    input:
        annotation1=expand("{ref_path}/{species}_{build}_{release}/annotation.gtf",
            species=species_list[0],
            build=build_list[0],
            release=release_list[0],
            ref_path=config['META']['reference-directory']),
        annotation2=expand("{ref_path}/{species}_{build}_{release}/annotation.gtf",
            species=species_list[1],
            build=build_list[1],
            release=release_list[1],
            ref_path=config['META']['reference-directory']),
    output:
        "{}/{}_{}_{}/annotation.gtf".format(
            config['META']['reference-directory'],
            species,
            build,
            release)
    params:
        build_list=build_list,
        release_list=release_list,
        species_list=species_list
    run:
        import datetime
        import re
        header1="#!Mixed reference of {} and {}\n".format(
            species_list[0],
            species_list[1])
        header2="#!genome-builds GRC{}{} GRC{}{}\n".format(
            species_list[0].lower()[0],
            build_list[0],
            species_list[1].lower()[0],
            build_list[1])
        header3="#!genome-releases {} {}\n".format(
            release_list[0],
            release_list[1])
        header4="#!genome-date {}\n".format(str(datetime.date.today()))
        header=[header1,header2,header3,header4]
        with open(input.annotation1[0]) as annotation1:
            with open(input.annotation2[0]) as annotation2:
                with open(output[0], 'w') as outfile:
                    outfile.writelines(header)
                    for line in annotation1:
                        if(not line.startswith('#!')):
                            outfile.write(re.sub('^',species_list[0],line))
                    for line in annotation2:
                        if(not line.startswith('#!')):
                            outfile.write(re.sub('^',species_list[1],line))


