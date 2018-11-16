import pysam
import re
import csv
from Bio import SeqIO
import gzip
from collections import defaultdict
import sys

#This function fills in a dict with readids
#and their corresponding cell and umi barcodes until it finds
#a specific read id

discard_secondary_alignements = snakemake.params['discard_secondary_alignements']

barcodes_struct = {
	'BC_start':snakemake.params['BC_start'],
	'BC_end':snakemake.params['BC_end'],
	'UMI_start':snakemake.params['UMI_start'],
	'UMI_end':snakemake.params['UMI_end']
	}


def parse_barcodes(fastq_parser, query_name, read_barcodes, barcodes_struct):
	for fastq_R1 in fastq_parser:
		read_barcodes[fastq_R1.id]['XC'] = str(fastq_R1.seq)[barcodes_struct['BC_start']:barcodes_struct['BC_end']]
		read_barcodes[fastq_R1.id]['XM'] = str(fastq_R1.seq)[barcodes_struct['UMI_start']:barcodes_struct['UMI_end']]
		if(read_barcodes[fastq_R1.id]['XM']==''):
			sys.SystemExit('UMI empty for read {}.\n The barcode is: {}.\nWhole entry is:{}'.format(fastq_R1.id, fastq_R1.seq,fastq_R1))
		if (fastq_R1.id == query_name):
			return(fastq_parser,read_barcodes)
	return(fastq_parser,read_barcodes)
	
infile_bam = pysam.AlignmentFile(snakemake.input[0], "rb")

fastq_parser = SeqIO.parse(gzip.open(snakemake.input[1], "rt"), "fastq")

outfile = pysam.AlignmentFile(snakemake.output[0], "wb", template=infile_bam)

read_barcodes = defaultdict(lambda :{'XC':'','XM':''})

for bam_read in infile_bam:
	if(discard_secondary_alignements & bam_read.is_secondary):
		continue
	if (bam_read.query_name) in read_barcodes:
		current_barcodes = read_barcodes.pop(bam_read.query_name)
		bam_read.set_tags([
			('XC', current_barcodes['XC'],'Z'),
			('XM', current_barcodes['XM'],'Z')])
	else:
		fastq_parser,read_barcodes = parse_barcodes(fastq_parser, bam_read.query_name, read_barcodes, barcodes_struct)
		if (bam_read.query_name) not in read_barcodes:
			raise SystemExit('Read from mapped file is missing in reference fastq file!')
			os.remove(snakemake.output[0])
		current_barcodes = read_barcodes.pop(bam_read.query_name)
		bam_read.set_tags([
			('XC', current_barcodes['XC'],'Z'),
			('XM', current_barcodes['XM'],'Z')])
	outfile.write(bam_read)