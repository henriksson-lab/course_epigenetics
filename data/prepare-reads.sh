#!/bin/bash
# Teacher script: build the small read sets used in the alignment, QC and
# ATAC-seq units. Students do not need to run this if the prepared files are
# provided to them (see data/README.md).
#
#   data/prepare-reads.sh            # needs ~25 GB of scratch, samtools, python3, curl
#
# For each ENCODE dataset below it downloads the original FASTQ files and the
# ENCODE alignment (BAM), takes the names of all reads that ENCODE aligned to
# chr21 or chrM, and keeps exactly those reads from the original FASTQ. The
# result is untrimmed, real sequencing output (adapters and all) that aligns to
# a small chr21 + chrM reference in minutes. The big downloads are deleted.
#
# Datasets (GM12878, GRCh38):
#   ATAC-seq  ENCSR095QNB (Snyder lab), paired-end 2x100, replicates 1 and 2
#   CTCF ChIP ENCSR000DZN (Snyder lab), single-end 36 bp, replicate 1
#   input     ENCSR000EYX, the control for the CTCF experiment, single-end 28 bp
set -euo pipefail
cd "$(dirname "$0")"
TMP=${TMP:-prepare-tmp}
mkdir -p "$TMP"
enc() { echo "https://www.encodeproject.org/files/$1/@@download/$1.$2"; }
get() { [ -s "$TMP/$2" ] || curl -fL -sS --retry 3 -o "$TMP/$2" "$1"; }

# names of reads aligned (primary) to chr21 or chrM in a BAM
names() {
  local bam=$TMP/$1
  samtools index -@ 8 "$bam" 2>/dev/null || { samtools sort -@ 8 -o "$bam.s" "$bam" && mv "$bam.s" "$bam" && samtools index "$bam"; }
  samtools view -@ 8 -F 0x900 "$bam" chr21 chrM | cut -f1 | sed 's,/[12]$,,' | sort -u -S 4G > "$TMP/$2"
}

# reference: chr21 + chrM (UCSC hg38)
if [ ! -s hg38_chr21_chrM.fa ]; then
  curl -fsS https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz | gunzip -c >  hg38_chr21_chrM.fa
  curl -fsS https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz  | gunzip -c >> hg38_chr21_chrM.fa
fi

# ATAC-seq, paired-end: rep  BAM(unfiltered)  R1  R2
while read rep bam r1 r2; do
  out=gm12878_atac_${rep}
  [ -s ${out}_R2.fastq.gz ] && { echo "have $out"; continue; }
  get "$(enc $bam bam)" $bam.bam; get "$(enc $r1 fastq.gz)" $r1.fastq.gz; get "$(enc $r2 fastq.gz)" $r2.fastq.gz
  names $bam.bam $out.names
  python3 -I filter_fastq_by_names.py "$TMP/$out.names" "$TMP/$r1.fastq.gz" ${out}_R1.fastq.gz
  python3 -I filter_fastq_by_names.py "$TMP/$out.names" "$TMP/$r2.fastq.gz" ${out}_R2.fastq.gz
  rm -f "$TMP/$bam.bam"* "$TMP/$r1.fastq.gz" "$TMP/$r2.fastq.gz"
done <<'EOF'
rep1 ENCFF290DWH ENCFF411MVZ ENCFF736ZWY
rep2 ENCFF878MKE ENCFF196VAZ ENCFF195FRB
EOF

# ChIP-seq and input, single-end: name  BAM(unfiltered)  FASTQ
while read out bam r1; do
  [ -s $out.fastq.gz ] && { echo "have $out"; continue; }
  get "$(enc $bam bam)" $bam.bam; get "$(enc $r1 fastq.gz)" $r1.fastq.gz
  names $bam.bam $out.names
  python3 -I filter_fastq_by_names.py "$TMP/$out.names" "$TMP/$r1.fastq.gz" $out.fastq.gz
  rm -f "$TMP/$bam.bam"* "$TMP/$r1.fastq.gz"
done <<'EOF'
gm12878_ctcf_chip ENCFF764YXZ ENCFF000VUW
gm12878_ctcf_input ENCFF906AMN ENCFF000VWV
EOF

ls -la gm12878_*fastq.gz hg38_chr21_chrM.fa
