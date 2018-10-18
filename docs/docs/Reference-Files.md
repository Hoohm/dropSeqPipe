Reference files generation
-----------------
Before running the pipeline you will need to download a reference genome as well as the GTF annotation.
dropSeqPipe is based on multiple files derived from those two files. Those will be automatically generated when you run the pipeline.

As an example, we use ensembl reference and annotation [located here](http://www.ensembl.org/info/data/ftp/index.html/). The fasta reference is the DNA and the GTF is the gene sets.

All you need to do is put the reference genome as well as the GTF file in a folder (extentions are crucial. it won't run otherwise):


```
genome.fasta
annotation.gtf
```

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