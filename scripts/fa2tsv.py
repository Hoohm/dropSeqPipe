#' ---
#' title:  fa2tsv.py
#' author: Sebastian Mueller (sebm_at_posteo.de)
#' date:   2019-03-04-
#' Convertes fasta file into tab seperated file suitable as input for FastQC.
#' This is to have FastQC using customized adapters using the -a option
#' ---

import sys
from Bio import SeqIO

number_bp = 12

with open( snakemake.output['tsv'], "w" ) as output:
    for seq_record in SeqIO.parse(snakemake.input['fa'], "fasta"):
        myline = (str(seq_record.id)) + "\t" + str(seq_record.seq[0:(number_bp + 1)]) + "\n"
        output.write(myline)

