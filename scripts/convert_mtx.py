#Converts the long format given by the dropseqtools v2.0.0 into the sparse mtx format.
#Output provides features, cell barcodes and counts in seperate files.
# Can handle one or multiple samples at a time

import os
import subprocess


samples = snakemake.params['samples']
barcodes = {}
features = {}

out_folder = os.path.dirname(snakemake.output['mtx'])
out_barcodes = snakemake.output.barcodes
out_features = snakemake.output.features
mtx = snakemake.output.mtx
temp_mtx = os.path.join(out_folder,'temp_umi.mtx')
header = os.path.join(out_folder,'header.mtx')
n_lines=0
barcode_index = 1
feature_index = 1

with open(temp_mtx,'w') as mtx_stream:
	for i,sample in enumerate(snakemake.input):
		if samples[i] not in sample:
			sys.exit("Sample name not found in file path")
		with open(sample,'r') as input_file:
			next(input_file) # skip first line
			for line in input_file:
				barcode,feature,count = line.strip().split('\t')
				if(not isinstance(samples,str)):
					barcode = samples[i] + '_' + barcode
				if(barcode not in barcodes):
					barcodes[barcode] = barcode_index
					barcode_index += 1
				if(feature not in features):
					features[feature] = feature_index
					feature_index += 1
				mtx_stream.write('{} {} {}\n'.format(features[feature],barcodes[barcode],count))
				n_lines +=1

with open(out_barcodes,'w') as barcodes_outfile:
	for barcode in barcodes:
		barcodes_outfile.write('{}\n'.format(barcode))

with open(out_features,'w') as features_outfile:
	for feature in features:
		features_outfile.write('{}\n'.format(feature))

with open(header,'w') as header_outfile:
	header_outfile.write("%%MatrixMarket matrix coordinate real general\n")
	header_outfile.write('{} {} {}\n'.format(len(features), len(barcodes), n_lines))

subprocess.call("cat {} {} > {}".format(header, temp_mtx, mtx), shell=True)

os.remove(temp_mtx)
os.remove(header)