#!/bin/bash
# The command-line part of the methylation-bisulfite exercise in one script:
# trim, build the bisulfite index for chr21, align with Bismark, deduplicate,
# extract methylation calls, merge the strands of each CpG.
#
#   bash run-bismark.sh OUTDIR        # about 5 minutes with 8 threads
#
# Needs fastp, Bismark (with bowtie2 and samtools) on the PATH, and the course
# data in $D (default: ../data relative to this script).
set -euo pipefail
OUT=${1:-bismark-run}
D=${D:-$(cd "$(dirname "$0")/../data" && pwd)}
mkdir -p "$OUT" && cd "$OUT"

fastp -i $D/gm12878_wgbs_reads.fastq.gz -o trimmed.fastq.gz \
  --adapter_sequence AGATCGGAAGAGC --length_required 30 -w 4 -j fastp.json -h fastp.html
mkdir -p genome && zcat $D/hg38_chr21.fa.gz > genome/chr21.fa
bismark_genome_preparation --parallel 4 genome
bismark --genome genome --parallel 2 -p 2 trimmed.fastq.gz
deduplicate_bismark --single --bam trimmed_bismark_bt2.bam
bismark_methylation_extractor --single-end --gzip --bedGraph --parallel 4 trimmed_bismark_bt2.deduplicated.bam
coverage2cytosine --merge_CpG --zero_based --genome_folder genome -o gm_cpg \
  trimmed_bismark_bt2.deduplicated.bismark.cov.gz
