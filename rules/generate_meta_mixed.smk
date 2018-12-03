rule rename_genome:
    input:
        temp("{ref_path}/{species}_{build}_{release}/original_genome.fa")
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

rule merge_annotation:
    input:
        annotation1=expand("{ref_path}/{species}_{build}_{release}/original_annotation.gtf",
            species=species_list[0],
            build=build_list[0],
            release=release_list[0],
            ref_path=config['META']['reference-directory']),
        annotation2=expand("{ref_path}/{species}_{build}_{release}/original_annotation.gtf",
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



#!genome-build GRCm38.p6
#!genome-version GRCm38
#!genome-date 2012-01
#!genome-build-accession NCBI:GCA_000001635.8
#!genebuild-last-updated 2018-06
