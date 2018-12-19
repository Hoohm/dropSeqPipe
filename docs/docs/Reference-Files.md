Reference files
-----------------
From version 0.4 on, reference files are automatically downloaded by the pipeline. Mixed references are also downloaded and merged automatically. Since sometimes you still want to use your own reference you can bypass the download by creating your own `genome.fa` and `annotation.gtf` file.

Snakemake generates file based on paths. If you want to use a custom reference you have to name it properly for snakemake to find it.

Here is an example:

Let's assume this is you configuration for the META section:
```
META:
    species:
        funky_species_name:
            build: A
            release: 1
    ratio: 0.2
    reference-directory: /absolute/path/to/references
    gtf_biotypes: gtf_biotypes.yaml
```

You need to provide the following files

```
/absolute/path/to/references/funky_species_name_A_1/genome.fa
/absolute/path/to/references/funky_species_name_A_1/annotation.gtf
```

This will stop dropSeqPipe from downloading a new reference.


Once the pipeline has run completely, the folder will look like this:

```
genome.fa
annotation.gtf
annotation.refFlat
annotation_reduced.gtf
genome.consensus_introns.intervals
genome.dict
genome.exons.intervals
genome.genes.intervals
genome.intergenic.intervals
genome.rRNA.intervals
STAR_INDEX/SA_read_length/
```

Note: The STAR index will be built based on the read length of your mRNA read (Read2).
If you have different lengths, it will produce multiple indexes.

Finally, you can now [run the pipeline](https://github.com/Hoohm/dropSeqPipe/wiki/Running-dropSeqPipe)