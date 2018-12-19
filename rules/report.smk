import os
import glob

localrules: create_publication_text

def get_yamls(wildcards):
    files = glob.glob('.snakemake/conda/*.yaml')
    return(files)

rule create_publication_text:
    input:
        config_file=configfile_path,
        yaml_files=get_yamls
    output:
        '{results_dir}/reports/publication_text.html'
    script:
        "../scripts/publication_text.Rmd"