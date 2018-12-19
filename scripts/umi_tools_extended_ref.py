import os
import sys
from collections import defaultdict
import pickle


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


mapping=defaultdict(dict)
barcode_ref = set()
barcode_ext_ref = set()


with open(snakemake.input['whitelist'],'r') as whitelist:
	for line in whitelist:
		if len(line.strip().split()) == 2:  # This means we didn't find any other linked barcode
			(reference,counts_ref) = line.strip().split()
			mapping[0][reference]= defaultdict()
			mapping[0][reference]['ref'] = reference
			mapping[0][reference]['count'] = 0
			mapping[0][reference]['lanes'] = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}
			barcode_ref.add(reference)
			continue
		(reference,extended_ref,counts_ref,counts_ext) = line.strip().split()
		mapping[0][reference]= defaultdict()
		mapping[0][reference]['ref'] = reference
		mapping[0][reference]['count'] = 0
		mapping[0][reference]['lanes'] = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}
		barcode_ref.add(reference)
		for barcode in extended_ref.split(','):
			mapping[1][barcode] = defaultdict()
			mapping[1][barcode]['ref'] = reference
			mapping[1][barcode]['count'] = 0
			mapping[1][barcode]['lanes'] = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}
		barcode_ext_ref.update(extended_ref.split(','))

# Save mapping and references to reuse later.
save_obj(obj=mapping, name=snakemake.output.barcode_mapping)
save_obj(obj=barcode_ref,name=snakemake.output.barcode_ref)
save_obj(obj=barcode_ext_ref,name=snakemake.output.barcode_ext_ref)