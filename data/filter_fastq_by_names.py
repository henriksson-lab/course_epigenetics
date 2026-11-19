#!/usr/bin/env python3
"""Keep the FASTQ records whose read name is in a list.

    python3 filter_fastq_by_names.py names.txt in.fastq.gz out.fastq.gz

The name is the first word of the header without the leading '@' and without
a trailing /1 or /2, which is how aligners write it into a BAM file.
"""
import gzip
import sys

names_file, fq_in, fq_out = sys.argv[1:4]
with open(names_file) as f:
    keep = {line.strip() for line in f if line.strip()}

n_in = n_out = 0
with gzip.open(fq_in, "rt") as fin, gzip.open(fq_out, "wt", compresslevel=6) as fout:
    while True:
        rec = [fin.readline() for _ in range(4)]
        if not rec[0]:
            break
        n_in += 1
        name = rec[0][1:].split()[0]
        if name.endswith("/1") or name.endswith("/2"):
            name = name[:-2]
        if name in keep:
            fout.writelines(rec)
            n_out += 1
print(f"{fq_in}: kept {n_out} of {n_in} reads", file=sys.stderr)
