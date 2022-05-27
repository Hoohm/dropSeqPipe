from itertools import combinations, product
from collections import defaultdict
from copy import deepcopy
import pickle
from shutil import copyfile


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def generate_all(barcode, reference, mapping, edit_distance):
    mutants = generate_mutants(barcode, edit_distance)
    for mutant in mutants:
        if(mutant not in reference):
            reference.add(mutant)
            mapping[edit_distance][mutant]['ref'] = barcode
            mapping[edit_distance][mutant]['count'] = 0
            mapping[edit_distance][mutant]['lanes'] = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}

    mapping['unknown']=defaultdict()
    return(reference, mapping)


def generate_mutants(sequence, d=1):
    """Taken from stackoverflow: https://stackoverflow.com/a/19823295/9178565"""
    N = len(sequence)
    letters = 'ACGTN'
    pool = list(sequence)
    for indices in combinations(range(N), d):
        for replacements in product(letters, repeat=d):
            skip = False
            for i, a in zip(indices, replacements):
                if pool[i] == a: skip = True
            if skip: continue

            keys = dict(zip(indices, replacements))
            yield ''.join([pool[i] if i not in indices else keys[i] 
                           for i in range(N)])

# Create empty sets and defaultdicts
barcode_ref = set()
mapping=defaultdict(dict)

# Initiate ref and mapping with the given barcodes
with open(snakemake.input.whitelist,'r') as ref_file:
    for line in ref_file.readlines():
        barcode = line.strip()
        barcode_ref.add(barcode)
        mapping[0][barcode]=defaultdict(dict)
        mapping[0][barcode]['ref'] = barcode
        mapping[0][barcode]['count'] = 0
        mapping[0][barcode]['lanes'] = {'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}

barcode_ext_ref = deepcopy(barcode_ref)
# For now edit distance is one, but can be extended to a higher number later on.
max_edit_distance = 1
for edit_distance in range(1,max_edit_distance+1):
    mapping[edit_distance]=defaultdict(dict)
    for barcode in mapping[0]:
        (barcode_ext_ref,mapping) = generate_all(barcode, barcode_ext_ref, mapping, edit_distance)

# Delete given barcodes out of new reference. This helps later on when running "repair_barcodes.py"
barcode_ref = set(mapping[0])
barcode_ext_ref.difference_update(barcode_ref)

# Save mapping and references to reuse later.
save_obj(obj=mapping, name=snakemake.output.barcode_mapping)
save_obj(obj=barcode_ref,name=snakemake.output.barcode_ref)
save_obj(obj=barcode_ext_ref,name=snakemake.output.barcode_ext_ref)


copyfile(snakemake.input['whitelist'], snakemake.output['barcodes'])