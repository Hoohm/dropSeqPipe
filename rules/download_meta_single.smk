from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

localrules:
	download_annotation,
	download_genome

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