#!/bin/bash
# Teacher script: build the small data sets for the ChIP-seq units
# (ChIP-seq: the assay; Peak calling; Annotation and differential binding).
# Students do not need to run this if the prepared files are provided to them.
#
#   data/prepare-chip.sh          # needs ~10 GB of scratch, samtools >= 1.10, python3, curl
#   TMP=/scratch/me data/prepare-chip.sh
#
# Part 1, reads (same recipe as prepare-reads.sh): download the original FASTQ
# and ENCODE's unfiltered GRCh38 alignment, take the names of the reads ENCODE
# aligned to chr21 or chrM, and keep exactly those reads from the FASTQ.
#   gm12878_ctcf_chip_rep2.fastq.gz   CTCF ChIP-seq ENCSR000DZN replicate 2 (Snyder lab, SE 36 bp, Phred+64)
#   gm12878_h3k36me3_chip.fastq.gz    H3K36me3 ChIP-seq ENCSR000AKE replicate 1 (Bernstein lab, SE 51 bp)
#   gm12878_h3k36me3_chip_rep2.fastq.gz  replicate 2 of the same experiment (SE 51 bp)
#   gm12878_h3k36me3_input.fastq.gz   its input ENCSR000AKJ replicate 1 (Bernstein lab, SE 51 bp)
#
# Part 2, filtered alignments for differential binding: ENCODE's filtered GRCh38
# BAMs (MAPQ >= 30, duplicates removed; ENCODE pipeline 1.5.1), cut to chr21.
#   h3k27ac_{gm12878,k562}_rep{1,2}.chr21.bam(.bai)   H3K27ac ENCSR000AKC (GM12878, 51 bp) and ENCSR000AKP (K562, 36 bp)
#   input_{gm12878,k562}.chr21.bam(.bai)              inputs ENCSR000AKJ rep 1 and ENCSR000AKY rep 1
#   h3k27ac_idxstats.tsv                              reads per chromosome in the whole filtered BAMs
# The big downloads are deleted as soon as they have been used.
set -euo pipefail
cd "$(dirname "$0")"
TMP=${TMP:-prepare-tmp}
mkdir -p "$TMP"
enc() { echo "https://www.encodeproject.org/files/$1/@@download/$1.$2"; }
get() { [ -s "$TMP/$2" ] || curl -fL -sS --retry 3 -o "$TMP/$2" "$1"; }

names() {  # names of reads aligned (primary) to chr21 or chrM
  local bam=$TMP/$1
  samtools index -@ 8 "$bam" 2>/dev/null || { samtools sort -@ 8 -o "$bam.s" "$bam" && mv "$bam.s" "$bam" && samtools index "$bam"; }
  samtools view -@ 8 -F 0x900 "$bam" chr21 chrM | cut -f1 | sed 's,/[12]$,,' | sort -u -S 4G > "$TMP/$2"
}

# ---- part 1: single-end reads ------------------------------------------------
while read out bam r1; do
  [ -s $out.fastq.gz ] && { echo "have $out"; continue; }
  get "$(enc $bam bam)" $bam.bam; get "$(enc $r1 fastq.gz)" $r1.fastq.gz
  names $bam.bam $out.names
  python3 -I filter_fastq_by_names.py "$TMP/$out.names" "$TMP/$r1.fastq.gz" $out.fastq.gz
  rm -f "$TMP/$bam.bam"* "$TMP/$r1.fastq.gz" "$TMP/$out.names"
done <<'EOF'
gm12878_ctcf_chip_rep2 ENCFF340XBU ENCFF000VUU
gm12878_h3k36me3_chip ENCFF378CVF ENCFF000ATN
gm12878_h3k36me3_chip_rep2 ENCFF508IAX ENCFF000ATJ
gm12878_h3k36me3_input ENCFF383MXY ENCFF000ARK
EOF

# ---- part 2: filtered BAMs, chr21 only ---------------------------------------
[ -s h3k27ac_idxstats.tsv ] || printf "sample\tchrom\tlength\tmapped\tunmapped\n" > h3k27ac_idxstats.tsv
while read out bam; do
  [ -s $out.chr21.bam.bai ] && { echo "have $out"; continue; }
  get "$(enc $bam bam)" $bam.bam
  samtools index -@ 8 "$TMP/$bam.bam"
  samtools idxstats "$TMP/$bam.bam" | awk -v s=$out -v OFS='\t' '{ print s, $0 }' >> h3k27ac_idxstats.tsv
  samtools view -@ 8 -b -o $out.chr21.bam "$TMP/$bam.bam" chr21
  samtools index $out.chr21.bam
  rm -f "$TMP/$bam.bam"*
done <<'EOF'
h3k27ac_gm12878_rep1 ENCFF269GKF
h3k27ac_gm12878_rep2 ENCFF201OHW
h3k27ac_k562_rep1 ENCFF121RHF
h3k27ac_k562_rep2 ENCFF907MNY
input_gm12878 ENCFF346JZP
input_k562 ENCFF226FKB
EOF

ls -la gm12878_ctcf_chip_rep2.fastq.gz gm12878_h3k36me3_*.fastq.gz h3k27ac_*.chr21.bam input_*.chr21.bam h3k27ac_idxstats.tsv
