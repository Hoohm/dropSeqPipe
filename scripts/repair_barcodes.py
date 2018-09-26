import pickle
import pysam

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


infile_bam = pysam.AlignmentFile(snakemake.input.bam, "rb")
outfile = pysam.AlignmentFile(snakemake.output.bam, "wb", template=infile_bam)

mapping = load_obj(snakemake.input.barcode_mapping)
barcode_ref = load_obj(snakemake.input.barcode_ref)
barcode_ext_ref = load_obj(snakemake.input.barcode_ext_ref)


for bam_read in infile_bam:
	barcode = bam_read.get_tag('XC')
	if barcode in barcode_ref:
		mapping[0][barcode]['count'] += 1
		outfile.write(bam_read)
		continue
	elif barcode in barcode_ext_ref:
		# The barcode is in our extended reference. Change the barcode to the original one
		reference_barcode = mapping[1][barcode]['ref']
		mapping[1][barcode]['count'] += 1
		bam_read.set_tags([('XC',reference_barcode,'Z')])
		outfile.write(bam_read)
		continue
	else:
		# If the barcode is not found in the extended ref, then don't modify it.
		outfile.write(bam_read)

save_obj(obj=mapping, name=snakemake.output.barcode_mapping_counts)
