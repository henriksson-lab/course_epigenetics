#!/bin/bash
# Teacher script: build the small data sets for the units "CUT&RUN and CUT&Tag"
# and "Combining marks into chromatin states". Students do not need to run it
# if the prepared files are provided to them (see data/README.md).
#
#   data/prepare-cuttag.sh        # needs ~15 GB of scratch, bowtie2, samtools >= 1.16, python3, curl
#   TMP=/scratch/me data/prepare-cuttag.sh
#
# 1. CUT&Tag in K562 (Kaya-Okur et al. 2020, GEO GSE145187, paired-end 2 x 25 bp):
#    H3K27me3 (SRR12246717), H3K4me3 (SRR11074254), IgG (SRR11923224); the
#    replicate-1 samples of the Henikoff-lab CUT&Tag tutorial. All read pairs are
#    aligned once to the full GRCh38 genome (prebuilt Bowtie2 index) and to
#    E. coli K-12 MG1655 (U00096.3); we keep the original FASTQ records of the
#    pairs that aligned to chr21 or chrM, or concordantly to E. coli (the
#    carry-over spike-in). Output: k562_cuttag_{h3k27me3,h3k4me3,igg}_R{1,2}.fastq.gz
#    and ecoli_U00096.3.fa.
# 2. ENCODE ChIP-seq alignments (GRCh38, filtered "alignments" BAMs of the default
#    analysis) cut to chr21: K562 H3K27me3 / H3K4me3 / input, replicate 1 (for the
#    comparison with CUT&Tag), and GM12878 H3K4me3, H3K27ac, H3K4me1, H3K36me3,
#    H3K27me3, H3K9me3 and input, replicates 1 + 2 merged (for ChromHMM).
#    Output: *_chr21.bam + .bai.
set -euo pipefail
cd "$(dirname "$0")"
TMP=${TMP:-prepare-tmp}
T=${T:-8}
SAMTOOLS=${SAMTOOLS:-samtools}
mkdir -p "$TMP"
get() { if [ ! -s "$2" ]; then curl -fL -sS --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"; fi; }
enc() { echo "https://www.encodeproject.org/files/$1/@@download/$1.$2"; }

# ---- 1. CUT&Tag reads -------------------------------------------------------
get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U00096.3&rettype=fasta&retmode=text" ecoli_U00096.3.fa

need_ct=0
for s in h3k27me3 h3k4me3 igg; do [ -s k562_cuttag_${s}_R2.fastq.gz ] || need_ct=1; done
if [ $need_ct = 1 ]; then
  # prebuilt Bowtie2 index of GRCh38 without alt contigs (3.7 GB download, ~4 GB unpacked)
  if [ ! -s "$TMP/GRCh38_noalt_as/GRCh38_noalt_as.1.bt2" ]; then
    get https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip "$TMP/GRCh38_noalt_as.zip"
    (cd "$TMP" && unzip -q -o GRCh38_noalt_as.zip && rm GRCh38_noalt_as.zip)
  fi
  [ -s "$TMP/ecoli.1.bt2" ] || bowtie2-build --threads 4 ecoli_U00096.3.fa "$TMP/ecoli" > /dev/null
fi

while read s run; do
  out=k562_cuttag_$s
  [ -s ${out}_R2.fastq.gz ] && { echo "have $out"; continue; }
  ena=ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}/0${run: -2}/$run
  get https://$ena/${run}_1.fastq.gz "$TMP/${run}_1.fastq.gz"
  get https://$ena/${run}_2.fastq.gz "$TMP/${run}_2.fastq.gz"
  # pairs with a primary alignment of either mate on chr21 or chrM (GRCh38: names chr21, chrM)
  bowtie2 -p $T --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
      -x "$TMP/GRCh38_noalt_as/GRCh38_noalt_as" -1 "$TMP/${run}_1.fastq.gz" -2 "$TMP/${run}_2.fastq.gz" \
      2> "$TMP/$run.hg38.log" \
    | $SAMTOOLS view -F 0x904 - | awk '$3 == "chr21" || $3 == "chrM" { print $1 }' > "$TMP/$out.names"
  # pairs aligned concordantly to E. coli (Henikoff-lab spike-in settings)
  bowtie2 -p $T --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant \
      --phred33 -I 10 -X 700 -x "$TMP/ecoli" -1 "$TMP/${run}_1.fastq.gz" -2 "$TMP/${run}_2.fastq.gz" \
      2> "$TMP/$run.ecoli.log" \
    | $SAMTOOLS view -F 0x904 -f 2 - | cut -f1 >> "$TMP/$out.names"
  sort -u -S 2G -o "$TMP/$out.names" "$TMP/$out.names"
  python3 -I filter_fastq_by_names.py "$TMP/$out.names" "$TMP/${run}_1.fastq.gz" ${out}_R1.fastq.gz
  python3 -I filter_fastq_by_names.py "$TMP/$out.names" "$TMP/${run}_2.fastq.gz" ${out}_R2.fastq.gz
  rm -f "$TMP/${run}_1.fastq.gz" "$TMP/${run}_2.fastq.gz"
done <<'EOF'
h3k27me3 SRR12246717
h3k4me3 SRR11074254
igg SRR11923224
EOF

# ---- 2. ENCODE alignments cut to chr21 --------------------------------------
# name  accession(s): several replicates are merged after cutting to chr21
while read name accs; do
  out=${name}_chr21.bam
  [ -s $out.bai ] && { echo "have $out"; continue; }
  parts=()
  for acc in ${accs//,/ }; do
    get "$(enc $acc bam)" "$TMP/$acc.bam"
    $SAMTOOLS index -@ 4 "$TMP/$acc.bam"
    $SAMTOOLS view -@ 4 -b -o "$TMP/$acc.chr21.bam" "$TMP/$acc.bam" chr21
    rm -f "$TMP/$acc.bam" "$TMP/$acc.bam.bai"
    parts+=("$TMP/$acc.chr21.bam")
  done
  if [ ${#parts[@]} = 1 ]; then mv "${parts[0]}" $out; else $SAMTOOLS merge -f -o $out "${parts[@]}"; rm -f "${parts[@]}"; fi
  $SAMTOOLS index $out
done <<'EOF2'
k562_h3k27me3_chip ENCFF392ZKG
k562_h3k4me3_chip ENCFF564SVK
k562_input ENCFF226FKB
gm12878_h3k4me3_chip ENCFF843BWY,ENCFF126HII
gm12878_h3k27ac_chip ENCFF269GKF,ENCFF201OHW
gm12878_h3k4me1_chip ENCFF047NLO,ENCFF385FLM
gm12878_h3k36me3_chip ENCFF353YPB,ENCFF677MAG
gm12878_h3k27me3_chip ENCFF633BHN,ENCFF565DCK
gm12878_h3k9me3_chip ENCFF306ENK,ENCFF889UPU
gm12878_input ENCFF708SHO,ENCFF834GHC
EOF2

ls -la k562_cuttag_*fastq.gz ecoli_U00096.3.fa *_chr21.bam
