import pysam
import re
import csv


with open(snakemake.input[1], mode='r') as infile:
	reader = csv.reader(infile, delimiter='\t')
	next(reader, None)
	mapping = {rows[0]:rows[1] for rows in reader}

infile = pysam.AlignmentFile(snakemake.input[0], "rb")

outfile = pysam.AlignmentFile(snakemake.output[0], "wb", template=infile)

pattern = '|'.join(sorted(re.escape(k) for k in mapping))

for s in infile:
	old_tag = s.get_tag('XC')
	new_tag = re.sub(pattern, lambda m: mapping.get(m.group(0).upper()), s.get_tag('XC'), flags=re.IGNORECASE)
	s.set_tag('XC', new_tag)
	outfile.write(s)
