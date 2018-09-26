from Bio import SeqIO
import gzip
from collections import Counter

fastq_parser = SeqIO.parse(gzip.open(snakemake.input.R1, "rt"), "fastq")
sequences = []
n=0
for fastq_R1 in fastq_parser:
		sequences.append(str(fastq_R1.seq))
		n+=1
		if(n==10000000):
			break
def parse_barcodes(fastq_parser):
	counts={}
	ranges = range(5,len(sequences[0]))

	for cell_bc_length in ranges:
		counts[cell_bc_length] = list()
		for fastq_R1 in sequences:
			counts[cell_bc_length].append(fastq_R1[0:cell_bc_length])
	return(counts)
	
	

counts = parse_barcodes(fastq_parser)

with open(snakemake.output[0], "w") as outfile:
	outfile.write('bc_length,first_counts\n')
	for cell_bc_length in counts:
		outfile.write('{},{}\n'.format(cell_bc_length,str(Counter(counts[cell_bc_length]).most_common(100)[0][1])))
