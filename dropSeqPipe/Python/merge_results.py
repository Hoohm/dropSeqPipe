"""Merge results."""
import yaml
import pandas as pd
import re
import os
import sys


def merge_read_counts(samples, extension, skip, tail, wdir):
    """Merge read counts."""
    readcounts_list = []
    for sample in samples:
        readcounts_list.append(pd.read_csv(
            "{}/summary/{}{}".format(wdir, sample, extension),
            skiprows=skip,
            skipfooter=tail,
            sep="\t", usecols=[0, 1],
            names=["GENE", sample],
            engine='python',
            index_col=0),)
    final = pd.concat(readcounts_list, axis=1)
    return(final)


def sorted_nicely(l):
    """Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key=alphanum_key)


def main():
    """Main."""
    wdir = sys.argv[1]
    configfile = os.path.join(wdir, 'config.yaml')
    with open(configfile) as yaml_file:
        yaml_config = yaml.load(yaml_file)
        samples = sorted_nicely(yaml_config['Samples'])# - yaml_config.get('Exclude', {'None':0}).keys())
    gene_name_matrix = merge_read_counts(samples, extension="", skip=0, tail=5, wdir=wdir)
    ensembl_matrix = merge_read_counts(samples, extension=".ReadsPerGene.out.tab", skip=4, tail=0, wdir=wdir)
    with open(os.path.join(wdir, "summary", "gene_name_matrix.tsv"), "w") as csv_file:
        gene_name_matrix.to_csv(csv_file, sep="\t")
    with open(os.path.join(wdir, "summary", "ensembl_matrix.tsv"), "w") as csv_file:
        ensembl_matrix.to_csv(csv_file, sep="\t")

main()
