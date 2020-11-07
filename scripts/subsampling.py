import pysam

from collections import namedtuple, Counter, defaultdict


bam_file = "/dss/dssfs02/lwp-dss-0001/pr62lo/pr62lo-dss-0000/di49tom/data/10kpbmcREAP/results/samples/pbmc_10k/final.bam"
cell_barcodes = "/dss/dssfs02/lwp-dss-0001/pr62lo/pr62lo-dss-0000/di49tom/data/10kpbmcREAP/results/samples/pbmc_10k/barcodes.csv"

strand = "+|-"

infile_bam = pysam.AlignmentFile(bam_file, "rb", threads=4)


barcodes = set()
with open(cell_barcodes, "r") as cell_barcodes:
    for line in cell_barcodes:
        barcodes.add(line.rstrip())

n_cells = len(barcodes)

sub_sampling_reads_per_cell_step = 10000

sub_sampling_reads_step = sub_sampling_reads_per_cell_step * n_cells


# categories

no_barcode = 0

unmapped = 0
mapped = 0
total_reads = 0

count_matrix = {}


for bam_read in infile_bam:
    if total_reads % sub_sampling_reads_step == 0:
        # write results to file
        print(
            "No barcodes:{}\nMapped:{}\nUnmapped:{}\nTotal reads:{}\nSummed reads:{}\n".format(
                no_barcode,
                mapped,
                unmapped,
                total_reads,
                no_barcode + mapped + unmapped,
            )
        )
        pass

    total_reads += 1
    if bam_read.is_unmapped:
        unmapped += 1
        continue
    tags = bam_read.get_tags()
    if "XC" in tags:
        if tags["XC"] not in barcodes:
            no_barcode += 1
            continue
        else:
            cell_barcode = tags["XC"]
            umi = tags["XM"]
            gene = tags["gn"]
            # If it's a cell
            if cell_barcode not in count_matrix:
                count_matrix[cell_barcode] = defaultdict(Counter)
            count_matrix[cell_barcode][gene][umi] += 1
            mapped += 1
