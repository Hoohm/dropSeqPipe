import csv

full_whitelist = set()
unrecognized_barcodes = set()

with open(snakemake.input.whitelist,'r') as whitelist:
    for line in whitelist:
        full_whitelist.add(line.strip())


with open(snakemake.input.top_barcodes, 'r') as top_barcodes, open(snakemake.output.filtered_barcodes, 'w') as filtered_barcodes:
    n = 0
    csv_iter = csv.reader(top_barcodes, delimiter='\t')
    next(csv_iter)
    for row in csv_iter:
        barcode = row[1].strip('-1')
        if barcode in full_whitelist:
            filtered_barcodes.write("{}\n".format(barcode))
            n+=1
        else:
            unrecognized_barcodes.add(barcode)
        if n%snakemake.params.num_cells==0:
            with open(snakemake.log.filter_log, 'w') as unknown_barcodes:
                for barcode in unrecognized_barcodes:
                    unknown_barcodes.write('{}\n'.format(barcode))
                print('Found {} known barcodes, exiting'.format(snakemake.params.num_cells))
                break
        barcode = row[1].strip('-1')

