#!/bin/bash
# Teacher script: build the Hi-C data used in the 3D-genome units
# (3D genome organisation; Hi-C: from reads to contact matrices;
#  Calling compartments, TADs and loops). Students do not need to run it if the
# prepared files are provided to them (see data/README.md).
#
#   data/prepare-hic.sh                  # all three parts
#   data/prepare-hic.sh matrices         # only part 1 (no aligner needed, ~10 min)
#
# Needs: python3 with cooler 0.10, hic-straw, pandas, numpy (see the conda line
# below); for part 2 also bwa, samtools, pairtools 1.1 and ~25 GB of temporary space.
#
#   conda create -n hic -c conda-forge -c bioconda python=3.11 cooler=0.10 \
#       cooltools=0.7 pairtools=1.1 bwa samtools pybigwig bioframe hic-straw
#
# 1. Matrices: GM12878 in situ Hi-C, MboI, Rao et al. 2014, all replicates combined
#    (ENCODE experiment ENCSR410MDC, reprocessed by ENCODE with Juicer on GRCh38,
#    file ENCFF555ISR, MAPQ >= 30, ~55 GB). We never download the whole .hic file:
#    hic-straw reads only the blocks it needs over HTTP. We take raw (unnormalised)
#    counts and write
#      gm12878_hic_chr19_chr21.mcool   chr19 and chr21, cis only, 5 kb .. 500 kb, ICE-balanced per chromosome
#      gm12878_hic_genome_1mb.cool     all autosomes + chrX, cis and trans, 1 Mb, ICE-balanced
#
# 2. Reads: the first 30 million read pairs of one sequencing run of the same
#    experiment (SRR1658570 = library HIC001, 2 x 101 bp), streamed from ENA.
#    They are aligned to the whole genome (bwa mem -SP5M) and parsed with
#    pairtools; the full-genome statistics are kept as
#      gm12878_hic_fullgenome_stats.txt   (pairtools dedup stats of all 30 M pairs)
#    and the original FASTQ records of the pairs whose two sides were placed on
#    chr21 (or were unmapped / multi-mapped, with the other side on chr21) are kept as
#      gm12878_hic_R1.fastq.gz, gm12878_hic_R2.fastq.gz
#    so that students can run the real pipeline against chr21 only, in minutes.
#
# Set BWA_INDEX to an existing bwa index of GRCh38 with UCSC names (chr1, ...)
# to skip building one (~1 h, 5 GB RAM, 5 GB disk). THREADS (default 8),
# NPAIRS (default 30000000) and TMP (default data/prepare-hic-tmp) can be set too.
set -euo pipefail
cd "$(dirname "$0")"
D=$(pwd)
TMP=${TMP:-$D/prepare-hic-tmp}
THREADS=${THREADS:-8}
NPAIRS=${NPAIRS:-30000000}
PART=${1:-all}
HIC_URL=https://encode-public.s3.amazonaws.com/2022/08/26/bf310b4b-bce7-41e0-91a2-fdec48d127fc/ENCFF555ISR.hic
mkdir -p "$TMP"

# ---- 1. matrices from the ENCODE .hic file ------------------------------------
cat > "$TMP/hic_to_cool.py" <<'EOF'
# Extract chromosomes from a remote .hic file into a .cool (raw counts) with hic-straw.
# usage: python hic_to_cool.py URL out.cool RESOLUTION chrA,chrB,... [--trans]
import sys
import numpy as np, pandas as pd, hicstraw, cooler
url, out, res, chroms = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4].split(",")
trans = "--trans" in sys.argv
h = hicstraw.HiCFile(url)
size = {c.name: c.length for c in h.getChromosomes()}
bins = cooler.binnify(pd.Series({c: size[c] for c in chroms}), res)
offset = {c: int(np.flatnonzero(bins.chrom.values == c)[0]) for c in chroms}
parts = []
for i, a in enumerate(chroms):
    for b in chroms[i:]:
        if a != b and not trans:
            continue
        r = h.getMatrixZoomData(a, b, "observed", "NONE", "BP", res).getRecords(0, size[a], 0, size[b])
        x = np.array([q.binX for q in r], dtype=np.int64) // res + offset[a]
        y = np.array([q.binY for q in r], dtype=np.int64) // res + offset[b]
        v = np.rint([q.counts for q in r]).astype(np.int32)
        parts.append(pd.DataFrame({"bin1_id": np.minimum(x, y), "bin2_id": np.maximum(x, y), "count": v}))
        print(a, b, len(r), int(v.sum()), file=sys.stderr)
pix = pd.concat(parts).groupby(["bin1_id", "bin2_id"], as_index=False)["count"].sum()
pix = pix[pix["count"] > 0].sort_values(["bin1_id", "bin2_id"])
cooler.create_cooler(out, bins, pix, ordered=True, assembly="hg38",
                     metadata={"source": url, "counts": "raw, MAPQ>=30, deduplicated (ENCODE Juicer)"})
EOF

if [ ! -s gm12878_hic_chr19_chr21.mcool ]; then
  echo "== chr19 + chr21 at 5 kb from $HIC_URL"
  python3 -I "$TMP/hic_to_cool.py" "$HIC_URL" "$TMP/chr19_chr21_5kb.cool" 5000 chr19,chr21
  cooler zoomify -p 4 -r 5000,10000,25000,50000,100000,250000,500000 \
    --balance --balance-args '--cis-only --nproc 4' -o gm12878_hic_chr19_chr21.mcool "$TMP/chr19_chr21_5kb.cool"
fi
if [ ! -s gm12878_hic_genome_1mb.cool ]; then
  echo "== genome-wide at 1 Mb"
  chroms=$(printf 'chr%s,' $(seq 1 22) X); chroms=${chroms%,}
  python3 -I "$TMP/hic_to_cool.py" "$HIC_URL" "$TMP/genome_1mb.cool" 1000000 "$chroms" --trans
  cooler balance "$TMP/genome_1mb.cool" && mv "$TMP/genome_1mb.cool" gm12878_hic_genome_1mb.cool
fi
[ "$PART" = matrices ] && exit 0

# ---- 2. reads -----------------------------------------------------------------
ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/000/SRR1658570/SRR1658570
for r in 1 2; do
  f=$TMP/SRR1658570_head_$r.fastq.gz
  [ -s "$f" ] || { echo "== stream the first $NPAIRS pairs of ${ENA##*/}_$r"
    # head stops reading early, so curl and zcat end with a broken pipe: ignore their exit status
    { curl -sS --retry 3 "${ENA}_$r.fastq.gz" 2>/dev/null || true; } | { zcat 2>/dev/null || true; } \
      | head -n $((4 * NPAIRS)) | gzip -1 > "$f.part" && mv "$f.part" "$f"; }
done

if [ -z "${BWA_INDEX:-}" ]; then
  BWA_INDEX=$TMP/hg38_main.fa
  if [ ! -s "$BWA_INDEX.bwt" ]; then
    echo "== build a bwa index of the primary hg38 chromosomes (about 1 h)"
    [ -s hg38.fa.bgz.fai ] || { echo "need data/hg38.fa.bgz (data/download.sh gene-regulation)"; exit 1; }
    samtools faidx hg38.fa.bgz $(printf 'chr%s ' $(seq 1 22) X Y M) > "$BWA_INDEX"
    bwa index "$BWA_INDEX"
  fi
fi
[ -s "$BWA_INDEX.fai" ] || samtools faidx "$BWA_INDEX"
cut -f1,2 "$BWA_INDEX.fai" > "$TMP/genome.chrom.sizes"

if [ ! -s "$TMP/full.pairs.gz" ]; then
  echo "== align to the whole genome and parse"
  bwa mem -SP5M -t "$THREADS" "$BWA_INDEX" \
      "$TMP/SRR1658570_head_1.fastq.gz" "$TMP/SRR1658570_head_2.fastq.gz" 2> "$TMP/bwa.log" \
    | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
        --drop-sam --drop-seq -c "$TMP/genome.chrom.sizes" -o "$TMP/full.pairs.gz"
fi
if [ ! -s gm12878_hic_fullgenome_stats.txt ]; then
  pairtools sort --nproc "$THREADS" --tmpdir "$TMP" "$TMP/full.pairs.gz" \
    | pairtools dedup --max-mismatch 3 --output-stats gm12878_hic_fullgenome_stats.txt -o /dev/null
fi

if [ ! -s gm12878_hic_R2.fastq.gz ]; then
  echo "== keep the pairs with both sides on chr21 (or one side on chr21, the other unplaced)"
  zcat "$TMP/full.pairs.gz" | grep -v '^#' \
    | awk '($2 == "chr21" || $2 == "!") && ($4 == "chr21" || $4 == "!") && ($2 == "chr21" || $4 == "chr21") { print $1 }' \
    > "$TMP/chr21.names"
  wc -l "$TMP/chr21.names"
  python3 -I filter_fastq_by_names.py "$TMP/chr21.names" "$TMP/SRR1658570_head_1.fastq.gz" gm12878_hic_R1.fastq.gz
  python3 -I filter_fastq_by_names.py "$TMP/chr21.names" "$TMP/SRR1658570_head_2.fastq.gz" gm12878_hic_R2.fastq.gz
fi

ls -la gm12878_hic_*
echo "done; remove $TMP when you are happy with the result"
