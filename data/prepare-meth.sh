#!/bin/bash
# Teacher script: build the small DNA-methylation data files used in
# *DNA methylation biology*, *Bisulfite and enzymatic sequencing* and
# *Differential methylation*. Students do not need to run this if the prepared
# files are provided to them (see data/README.md).
#
#   data/prepare-meth.sh        # needs ~3 GB temporary space, curl, python3,
#                               # samtools >= 1.10, bgzip/tabix, Bismark; ~25 min
#
# 1. Per-CpG methylation calls on chr21 (ENCODE WGBS, GRCh38, gemBS pipeline,
#    "methylation state at CpG", bed9+), two replicates each of GM12878 and K562.
#    The genome-wide files (~0.6 GB each) are downloaded, chr21 is kept.
# 2. Bisulfite reads for the alignment exercise: the read 1 sequences that
#    ENCODE aligned to chr21:33,000,001-35,000,000 in GM12878 WGBS replicate 1
#    (ENCSR890UQO, alignments ENCFF072CRF, from FASTQ ENCFF798RSS, 125 bp).
#    The 43 GB BAM is not downloaded: bam_slice_remote.py fetches only the
#    ~30 MB that hold the region (HTTP range requests). samtools fastq turns
#    the alignments back into the original reads (reverse-complementing reads
#    on the minus strand), unaligned and untrimmed. Reads are shuffled so that
#    `head` gives a random subsample.
# 3. A paired-end M-bias table for the lecture (read 1 and read 2): the first
#    800 MB of ENCODE's older Bismark alignment of the same replicate
#    (ENCFF681ASN, paired-end, unsorted, so a random sample of ~3 M pairs),
#    run through bismark_methylation_extractor --mbias_only. Needs Bismark
#    (set BISMARK=/path/to/Bismark-0.25.1 if it is not on the PATH); ~10 min.
set -euo pipefail
cd "$(dirname "$0")"
TMP=${TMP:-prepare-tmp}
mkdir -p "$TMP"
SAMTOOLS=${SAMTOOLS:-samtools}
enc() { echo "https://www.encodeproject.org/files/$1/@@download/$1.$2"; }

# 1. per-CpG calls, chr21 only
while read out acc; do
  [ -s $out ] && { echo "have $out"; continue; }
  echo "get $acc -> $out"
  [ -s "$TMP/$acc.bed.gz" ] || curl -fL -sS --retry 3 -o "$TMP/$acc.bed.gz" "$(enc $acc bed.gz)"
  zcat "$TMP/$acc.bed.gz" | awk -F'\t' '$1 == "chr21"' | bgzip > $out
  tabix -f -p bed $out
  rm -f "$TMP/$acc.bed.gz"
done <<'EOF'
gm12878_wgbs_rep1_chr21.bed.gz ENCFF570TIL
gm12878_wgbs_rep2_chr21.bed.gz ENCFF614QHA
k562_wgbs_rep1_chr21.bed.gz ENCFF660IHA
k562_wgbs_rep2_chr21.bed.gz ENCFF328NMN
EOF

# 2. bisulfite reads from one region of chr21
out=gm12878_wgbs_reads.fastq.gz
if [ ! -s $out ]; then
  url=$(curl -sI "$(enc ENCFF072CRF bam)" | grep -i '^location:' | tr -d '\r' | cut -d' ' -f2)
  python3 -I bam_slice_remote.py "$url" chr21 33000001 35000000 "$TMP/region.bam"
  $SAMTOOLS fastq -F 0x900 "$TMP/region.bam" 2> /dev/null |
    paste - - - - | shuf --random-source=<(yes) | tr '\t' '\n' | gzip > $out
  rm -f "$TMP/region.bam"
fi

# 3. paired-end M-bias table
out3=gm12878_wgbs_pe_mbias.txt
if [ ! -s $out3 ]; then
  url=$(curl -sI "$(enc ENCFF681ASN bam)" | grep -i '^location:' | tr -d '\r' | cut -d' ' -f2)
  for i in 0 1 2 3 4 5 6 7; do
    s=$((i * 100000000)); curl -fsS --retry 3 -r $s-$((s + 99999999)) -o "$TMP/part$i" "$url" &
  done; wait
  cat "$TMP"/part[0-7] > "$TMP/part.bam"; rm -f "$TMP"/part[0-7]
  # the file is cut off mid-record (samtools complains, hence || true): drop the last two records
  { $SAMTOOLS view -h "$TMP/part.bam" 2> /dev/null || true; } | head -n -2 | $SAMTOOLS view -b -o "$TMP/pe.bam" -
  ( cd "$TMP" && ${BISMARK:+$BISMARK/}bismark_methylation_extractor --paired-end --mbias_only --parallel 4 pe.bam > mbias.log 2>&1 )
  mv "$TMP/pe.M-bias.txt" $out3
  rm -f "$TMP"/part.bam "$TMP"/pe.bam "$TMP"/pe_splitting_report.txt
fi
ls -la gm12878_wgbs_*chr21.bed.gz k562_wgbs_*chr21.bed.gz $out $out3
