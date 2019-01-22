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
unknown_barcodes = set()

for bam_read in infile_bam:
	barcode = bam_read.get_tag('XC')
	#lane_number = bam_read.query_name.split(':')[3]
	if barcode in barcode_ref:
		mapping[0][barcode]['count'] += 1
		#mapping[0][barcode]['lanes'][lane_number] += 1
		outfile.write(bam_read)
		continue
	elif barcode in barcode_ext_ref:
		# The barcode is in our extended reference. Change the barcode to the original one
		reference_barcode = mapping[1][barcode]['ref']
		mapping[1][barcode]['count'] += 1
		#mapping[1][barcode]['lanes'][lane_number] += 1
		bam_read.set_tag('XC',reference_barcode,value_type='Z',replace=True)
		outfile.write(bam_read)
		continue
	else:
		# If the barcode is not found in the extended ref, then don't modify it.
		if barcode in unknown_barcodes:
			mapping['unknown'][barcode]['count'] += 1
			#mapping['unknown'][barcode]['lanes'][lane_number] += 1
		else:
			#mapping['unknown'][barcode] = {'count':1, 'lanes':{'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}}
			mapping['unknown'][barcode] = {'count':1}
			#mapping['unknown'][barcode]['lanes'][lane_number] += 1
			unknown_barcodes.add(barcode)
		outfile.write(bam_read)

save_obj(obj=mapping, name=snakemake.output.barcode_mapping_counts)
