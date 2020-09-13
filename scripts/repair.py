import gzip
import sys

from collections import namedtuple

from Bio import SeqIO, bgzf

# Used to convert the fastq stream into a file handle
from io import StringIO
from gzip import open as gzopen

print("started")


class Logs:
    def __init__(self):
        self.total_reads = 0
        self.R1_too_short = 0
        self.R2_too_short = 0
        self.both_too_short = 0
        self.written_reads = 0


logs = Logs()


with gzip.open(snakemake.input.R1, "rt") as r1_handle:
    with gzip.open(snakemake.input.R2, "rt") as r2_handle:
        with bgzf.BgzfWriter(snakemake.output.R2, "wb") as outgz_R2:
            with bgzf.BgzfWriter(snakemake.output.R1, "wb") as outgz_R1:
                parser_r1 = SeqIO.parse(r1_handle, "fastq")
                parser_r2 = SeqIO.parse(r2_handle, "fastq")
                for r1, r2 in zip(parser_r1, parser_r2):
                    r1_ok = False
                    r2_ok = False

                    if len(r2.seq) >= snakemake.params.min_r2_length:
                        r2_ok = True

                    if len(r1.seq) >= snakemake.params.barcode_length:
                        r1_ok = True

                    if r1_ok and r2_ok:
                        SeqIO.write(sequences=r2, handle=outgz_R2, format="fastq")
                        SeqIO.write(sequences=r1, handle=outgz_R1, format="fastq")
                        logs.written_reads += 1
                        logs.total_reads += 1
                        if logs.written_reads % 1000000 == 0:
                            print(
                                "Written {} reads from {} total reads".format(
                                    logs.written_reads, logs.total_reads
                                )
                            )
                    elif r1_ok and not r2_ok:
                        logs.total_reads += 1
                        logs.R2_too_short += 1
                    elif not r1_ok and r2_ok:
                        logs.total_reads += 1
                        logs.R1_too_short += 1
                    else:
                        logs.total_reads += 1
                        logs.both_too_short += 1

with open(snakemake.log.stats, "w") as logfile:
    logfile.write(
        "{},{},{},{},{}\n".format(
            "total_reads",
            "R1_toot_short",
            "R2_too_short",
            "both_too_short",
            "written_reads",
        )
    )
    logfile.write(
        "{},{},{},{},{}".format(
            logs.total_reads,
            logs.R1_too_short,
            logs.R2_too_short,
            logs.both_too_short,
            logs.written_reads,
        ),
    )
