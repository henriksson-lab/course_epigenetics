#!/bin/bash
# Teacher script: build the derived files of the integration block
# (*Linking regulatory elements to genes*, *Epigenomics and genetic variation*).
# Students do not need to run it if the files are provided.
#
#   data/prepare-integ.sh        # ~15 min, ~4 GB temporary space; needs curl, unzip,
#                                # bcftools, bedtools, python3, R with data.table,
#                                # GenomicRanges and rtracklayer
#
# Makes, in data/:
#   k562_abc_elements.tsv.gz       candidate regulatory elements in K562 (ATAC peaks
#                                  +-250 bp around summits, merged, plus protein-coding
#                                  promoters) with mean ATAC, DNase and H3K27ac signal
#   gwas_catalog_gws.tsv.gz        GWAS Catalog associations with p <= 5e-8 (one SNP,
#                                  GRCh38 position), trimmed to 10 columns
#   gwas_autoimmune_gm12878_alleles.tsv
#                                  reference/alternative alleles (dbSNP, via the NCBI
#                                  Variation API) of autoimmune GWAS SNPs inside GM12878
#                                  ATAC peaks, plus rs12740374 (SORT1)
#   na12878_chr21_het_snps.vcf.gz  heterozygous biallelic SNVs of GM12878 (= NA12878) on
#                                  chr21, GIAB v4.2.1 benchmark calls (+ .tbi index)
# Needs, already in data/ (data/download.sh): gm12878_atac_peaks.bed.gz,
# k562_atac_peaks.bed.gz, gencode.v46.basic.annotation.gtf.gz, hg38_blacklist.bed.gz.
set -euo pipefail
cd "$(dirname "$0")"
TMP=${TMP:-prepare-integ-tmp}
mkdir -p "$TMP"
enc() { echo "https://www.encodeproject.org/files/$1/@@download/$1.$2"; }
get() { [ -s "$TMP/$2" ] || curl -fL -sS --retry 3 -o "$TMP/$2" "$1"; }

# ---- 1. K562 candidate elements with ATAC, DNase and H3K27ac signal ---------------
# ATAC-seq ENCSR483RKN (Snyder; fold change, pooled reps), DNase-seq ENCSR000EOT
# (Stamatoyannopoulos; read-depth normalised, rep 1), H3K27ac ENCSR000AKP (Bernstein;
# fold change over input, reps 1+2). Same ENCODE pipelines as the GM12878 course data.
if [ ! -s k562_abc_elements.tsv.gz ]; then
  get "$(enc ENCFF754EAC bigWig)" k562_atac_fc.bigWig
  get "$(enc ENCFF414OGC bigWig)" k562_dnase_rdn.bigWig
  get "$(enc ENCFF849TDM bigWig)" k562_h3k27ac_fc.bigWig
  TMP="$TMP" Rscript - <<'EOF'
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(rtracklayer)})
tmp <- Sys.getenv("TMP"); chroms <- paste0("chr", c(1:22, "X"))
gtf <- import("gencode.v46.basic.annotation.gtf.gz", feature.type = "gene")
g <- gtf[gtf$gene_type == "protein_coding" & seqnames(gtf) %in% chroms]
prom <- resize(resize(g, 1, fix = "start"), 500, fix = "center")
p <- fread("k562_atac_peaks.bed.gz")[V1 %in% chroms]
pk <- GRanges(p$V1, IRanges(p$V2 + p$V10 + 1 - 250, width = 500))     # summit +- 250 bp
el <- reduce(c(granges(pk), granges(prom)), ignore.strand = TRUE)
el <- el[!overlapsAny(el, import("hg38_blacklist.bed.gz"))]
signal <- function(bw, gr) {                                            # mean signal per range
  cv <- import(bw, which = reduce(gr), as = "RleList"); s <- numeric(length(gr))
  for (ch in unique(as.character(seqnames(gr)))) { i <- which(as.character(seqnames(gr)) == ch)
    s[i] <- viewMeans(Views(cv[[ch]], ranges(gr)[i])) }
  s }
wide <- resize(el, width(el) + 500, fix = "center")                     # H3K27ac flanks the open site
out <- data.table(chr = as.character(seqnames(el)), start = start(el) - 1L, end = end(el),
                  promoter = overlapsAny(el, prom),
                  atac    = signal(file.path(tmp, "k562_atac_fc.bigWig"), el) * width(el) / 500,
                  dnase   = signal(file.path(tmp, "k562_dnase_rdn.bigWig"), el) * width(el) / 500,
                  h3k27ac = signal(file.path(tmp, "k562_h3k27ac_fc.bigWig"), wide) * width(wide) / 1000)
out[, (5:7) := lapply(.SD, round, 3), .SDcols = 5:7]
fwrite(out, "k562_abc_elements.tsv.gz", sep = "\t")
cat("K562 elements:", nrow(out), "\n")
EOF
  rm -f "$TMP"/k562_*.bigWig
fi

# ---- 2. GWAS Catalog: genome-wide significant associations -------------------------
if [ ! -s gwas_catalog_gws.tsv.gz ]; then
  get https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip gwas.zip
  unzip -p "$TMP/gwas.zip" | awk -F'\t' -v OFS='\t' '
    NR == 1 { print "snp", "chr", "pos", "risk_allele", "context", "trait", "trait_uri", "mlog10p", "study", "pubmed"; next }
    $28 != "" && $28 + 0 <= 5e-8 && $12 ~ /^([0-9]+|X)$/ && $13 ~ /^[0-9]+$/ && $22 ~ /^rs[0-9]+$/ {
      ra = $21; sub(/^[^-]*-/, "", ra)
      print $22, "chr" $12, $13, ra, $25, $35, $36, $29, $37, $2 }' | gzip > gwas_catalog_gws.tsv.gz
  rm -f "$TMP/gwas.zip"
fi

# ---- 3. alleles of autoimmune GWAS SNPs inside GM12878 ATAC peaks -------------------
if [ ! -s gwas_autoimmune_gm12878_alleles.tsv ]; then
  zcat gwas_catalog_gws.tsv.gz | awk -F'\t' -v OFS='\t' '
    $6 == "rheumatoid arthritis" || $6 == "systemic lupus erythematosus" || $6 == "multiple sclerosis" ||
    $6 == "inflammatory bowel disease" || $6 == "Crohn disease" || $6 == "ulcerative colitis" { print $2, $3 - 1, $3, $1 }' \
    | sort -k1,1 -k2,2n -u | bedtools intersect -u -a - -b gm12878_atac_peaks.bed.gz | cut -f4 | sort -u > "$TMP/ids.txt"
  echo rs12740374 >> "$TMP/ids.txt"
  python3 -I - "$TMP/ids.txt" > gwas_autoimmune_gm12878_alleles.tsv <<'EOF'
# GRCh38 position and dbSNP alleles of each rsID (NCBI Variation Services API, ~2 per second)
import json, sys, time, urllib.request
acc2chr = {f"NC_0000{i:02d}": f"chr{i}" for i in range(1, 23)}; acc2chr["NC_000023"] = "chrX"
print("snp\tchr\tpos\tref\talt")
for rs in (l.strip() for l in open(sys.argv[1])):
    if not rs: continue
    d = None
    for attempt in range(3):
        try:
            d = json.load(urllib.request.urlopen("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + rs[2:], timeout=60)); break
        except Exception:
            time.sleep(3)
    time.sleep(0.4)
    if not d or "primary_snapshot_data" not in d: continue
    for pl in d["primary_snapshot_data"]["placements_with_allele"]:
        if not pl.get("is_ptlp"): continue
        ref, alts, chrom, pos = None, [], None, None
        for al in pl["alleles"]:
            s = al["allele"].get("spdi")
            if not s: continue
            chrom = acc2chr.get(s["seq_id"].split(".")[0]); pos = s["position"] + 1
            if s["deleted_sequence"] == s["inserted_sequence"]: ref = s["deleted_sequence"]
            elif len(s["inserted_sequence"]) == 1: alts.append(s["inserted_sequence"])
        if chrom and ref and len(ref) == 1 and alts:
            print(f"{rs}\t{chrom}\t{pos}\t{ref}\t{','.join(alts)}")
EOF
  rm -f "$TMP/ids.txt"
fi

# ---- 4. GM12878 (NA12878) heterozygous SNVs on chr21, GIAB v4.2.1 -----------------
if [ ! -s na12878_chr21_het_snps.vcf.gz ]; then
  GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38
  bcftools view -r chr21 "$GIAB/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" -Ou \
    | bcftools view -g het -v snps -m2 -M2 -Ou | bcftools annotate -x INFO,^FORMAT/GT -Oz -o na12878_chr21_het_snps.vcf.gz
  bcftools index -t na12878_chr21_het_snps.vcf.gz
  rm -f HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
fi
rmdir "$TMP" 2>/dev/null || true
echo done
