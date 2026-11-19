#!/usr/bin/env bash
# Prepare the data for the units "Methylation arrays (Illumina EPIC)",
# "Epigenetic clocks and biomarkers" and "Methylation and chromatin from long reads".
#
#   data/prepare-meth2.sh              # everything (about 45 minutes, ~4 GB in data/, ~3 GB temporary)
#   data/prepare-meth2.sh arrays       # GSE42861 IDATs (151 samples) + sample table   (~1.3 GB)
#   data/prepare-meth2.sh clocks       # Hannum GSE40279 subset, Horvath coefficients, Moss atlas
#   data/prepare-meth2.sh longreads    # nanopore modBAM regions, trio SNPs, GM12878 WGBS subset
#
# Needs: curl, awk, samtools >= 1.15 (/usr/bin/samtools here), bcftools, tabix,
#        Rscript with data.table (clocks section), and data/hg38.fa.bgz
#        (built by data/download.sh; used to decode the nanopore CRAM files).
# Temporary files go to $TMPDIR (default /tmp).
set -euo pipefail
cd "$(dirname "$0")"
D=$(pwd)
TMP=$(mktemp -d -p "${TMPDIR:-/tmp}" meth2.XXXXXX)
trap 'rm -rf "$TMP"' EXIT
SAMTOOLS=${SAMTOOLS:-samtools}
what=${1:-all}

# ---------------------------------------------------------------------------
# 1. Methylation arrays: GSE42861 (Liu et al. 2013), 450k, peripheral blood
#    leukocytes of rheumatoid arthritis patients and controls. We take the 151
#    samples on the 14 BeadChips that carry both patients and controls.
# ---------------------------------------------------------------------------
arrays() {
  echo "== arrays: GSE42861 sample table and IDATs"
  curl -s https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/matrix/GSE42861_series_matrix.txt.gz \
    | gzip -dc 2>/dev/null | awk '/^!series_matrix_table_begin/ { exit } { print }' > "$TMP/header.txt" || true
  chips="7800246040 7766130166 7800246034 5730053006 5730053010 5730053011 5730053027 5730192057 7766130158 7800246024 7800246044 7800246057 7800246071 7800246188"
  awk -F'\t' -v chips="$chips" '
    function strip(s) { gsub(/"/, "", s); return s }
    BEGIN { n = split(chips, c, " "); for (i = 1; i <= n; i++) keep[c[i]] = 1 }
    $1 == "!Sample_geo_accession" { for (i = 2; i <= NF; i++) gsm[i] = strip($i); nf = NF }
    $1 == "!Sample_characteristics_ch1" {
      for (i = 2; i <= NF; i++) { v = strip($i); k = v; sub(/:.*/, "", k); sub(/^[^:]*: /, "", v); ch[k, i] = v } }
    $1 == "!Sample_supplementary_file" && !seen++ { for (i = 2; i <= NF; i++) url[i] = strip($i) }
    END {
      print "gsm\tbasename\tslide\tarray\tstatus\tage\tsex\tsmoking"
      for (i = 2; i <= nf; i++) {
        f = url[i]; sub(/.*\//, "", f); sub(/_Grn\.idat\.gz$/, "", f)          # GSM..._7800246024_R01C01
        split(f, p, "_")
        if (!(p[2] in keep)) continue
        st = (ch["subject", i] == "Patient") ? "RA" : "control"
        print gsm[i] "\t" f "\t" p[2] "\t" p[3] "\t" st "\t" ch["age", i] "\t" toupper(ch["gender", i]) "\t" ch["smoking status", i]
      } }' "$TMP/header.txt" > gse42861_samples.tsv
  echo "   $(($(wc -l < gse42861_samples.tsv) - 1)) samples in gse42861_samples.tsv"
  mkdir -p gse42861_idat
  tail -n +2 gse42861_samples.tsv | cut -f1,2 | while read -r gsm base; do
    for col in Grn Red; do
      f=gse42861_idat/${base}_${col}.idat.gz
      [ -s "$f" ] && continue
      curl -s --retry 5 -o "$f.part" \
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/${gsm}/suppl/${base}_${col}.idat.gz" && mv "$f.part" "$f"
    done
  done
  echo "   $(ls gse42861_idat/*.idat.gz | wc -l) IDAT files in gse42861_idat/"
  echo "== arrays: 450k probe annotation on GRCh38 and probe masks (Zhou, Laird & Shen 2017)"
  Z=https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450
  curl -sL -o hm450_hg38_manifest.tsv.gz $Z/HM450.hg38.manifest.tsv.gz
  curl -sL -o hm450_hg38_mask.tsv.gz $Z/HM450.hg38.mask.tsv.gz
}

# ---------------------------------------------------------------------------
# 2. Clocks: Hannum et al. 2013 (GSE40279, 656 whole-blood 450k arrays, ages
#    19-101) as a CpG subset; Horvath 2013 coefficients; Moss et al. 2018 atlas.
# ---------------------------------------------------------------------------
clocks() {
  echo "== clocks: Horvath coefficients and the Moss 2018 reference atlas"
  curl -sL -o horvath2013_coefficients.csv \
    "https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv"
  curl -sL -o moss2018_reference_atlas.csv \
    https://raw.githubusercontent.com/nloyfer/meth_atlas/master/reference_atlas.csv

  echo "== clocks: Hannum GSE40279 sample table"
  curl -s https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz \
    | gzip -dc 2>/dev/null | awk '/^!series_matrix_table_begin/ { exit } { print }' > "$TMP/hannum_header.txt" || true
  awk -F'\t' '
    function strip(s) { gsub(/"/, "", s); return s }
    $1 == "!Sample_geo_accession" { for (i = 2; i <= NF; i++) gsm[i] = strip($i); nf = NF }
    $1 == "!Sample_source_name_ch1" { for (i = 2; i <= NF; i++) src[i] = strip($i) }
    $1 == "!Sample_characteristics_ch1" {
      for (i = 2; i <= NF; i++) { v = strip($i); k = v; sub(/:.*/, "", k); sub(/^[^:]*: /, "", v); ch[k, i] = v } }
    END {
      print "gsm\tsample\tage\tsex\tethnicity\tsource\tplate"
      for (i = 2; i <= nf; i++)
        print gsm[i] "\t" src[i] "\t" ch["age (y)", i] "\t" ch["gender", i] "\t" ch["ethnicity", i] "\t" ch["source", i] "\t" ch["plate", i] }' \
    "$TMP/hannum_header.txt" > hannum_gse40279_samples.tsv
  echo "   $(($(wc -l < hannum_gse40279_samples.tsv) - 1)) samples"

  echo "== clocks: Hannum beta matrix (1.1 GB download) -> CpG subset"
  curl -s --retry 5 -o "$TMP/hannum_beta.txt.gz" \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/GSE40279_average_beta.txt.gz
  Rscript --vanilla - "$TMP/hannum_beta.txt.gz" hannum_gse40279_samples.tsv horvath2013_coefficients.csv \
      moss2018_reference_atlas.csv hannum_gse40279_beta_subset.tsv.gz <<'EOF'
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
b <- fread(a[1])                                           # ID_REF + one column per sample
setnames(b, 1, "cpg")
smp <- fread(a[2])
stopifnot(all(smp$sample %in% names(b)))
b <- b[grepl("^cg", cpg)]                                  # CpG probes only (drop ch. and rs probes)
m <- as.matrix(b[, smp$sample, with = FALSE]); rownames(m) <- b$cpg
cat("   full matrix:", nrow(m), "CpGs x", ncol(m), "samples\n")
# Keep the 20 000 most variable CpGs (age is NOT used to choose them), plus the
# CpGs of the Horvath clock and of the Moss atlas, so these can be applied, and
# two classic age CpGs from other studies (ELOVL2, FHL2) for plotting.
v <- apply(m, 1, var, na.rm = TRUE)
top <- names(sort(v, decreasing = TRUE))[1:20000]
hc <- fread(a[3], skip = 2)$CpGmarker
atl <- fread(a[4])$CpGs
lit <- c("cg16867657", "cg06639320")                       # ELOVL2 and FHL2, classic age CpGs
keep <- union(top, intersect(c(hc, atl, lit), rownames(m)))
out <- data.table(cpg = keep, round(m[keep, ], 4))
names(out) <- c("cpg", smp$gsm)                            # columns = GEO sample accessions
fwrite(out, a[5], sep = "\t", na = "NA")
cat("   kept", nrow(out), "CpGs:", length(top), "most variable,",
    sum(hc %in% keep), "of", length(hc) - 1, "Horvath CpGs,", sum(atl %in% keep), "atlas CpGs\n")
EOF
}

# ---------------------------------------------------------------------------
# 3. Long reads: HG001 (= NA12878, the donor of the GM12878 cell line) nanopore
#    R10.4.1 reads with 5mC/5hmC calls (ONT open data, giab_2023.05, 4 flow
#    cells), trio SNPs from 1000 Genomes, GM12878 WGBS (ENCODE) for comparison,
#    and EcoGII-treated (m6A) nanopore reads of HG002 (ONT chrom_acc_2025.06).
# ---------------------------------------------------------------------------
REGIONS="chr11:1960000-2050000 chr11:2690000-2710000 chr15:24940000-24970000 chrX:53150000-53340000 chr21:33000000-34000000"
regions_bed() { for r in $REGIONS; do echo "$r" | awk -F'[:-]' '{ print $1 "\t" $2 - 1 "\t" $3 }'; done; }

longreads() {
  [ -s hg38.fa.bgz.fai ] || { echo "need data/hg38.fa.bgz (+ .fai, .gzi): run data/download.sh first"; exit 1; }
  echo "== longreads: HG001 nanopore modBAM for $REGIONS"
  regions_bed > longread_regions.bed
  O=https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/hg001/sup
  parts=()
  if [ ! -s hg001_ont_regions.bam.bai ]; then
  for fc in PAO98424_a66982d8 PAO98424_ae91b84c PAO99309_94e07fab PAO99309_eece804b; do
    mkdir -p "$TMP/$fc"            # htslib caches remote indexes in the working directory
    ( cd "$TMP/$fc" && $SAMTOOLS view -b -M -L "$D/longread_regions.bed" -T "$D/hg38.fa.bgz" \
        -o "$TMP/$fc.bam" "$O/$fc.pass.cram" 2> "$TMP/$fc.log" ) || true
    $SAMTOOLS quickcheck "$TMP/$fc.bam" && parts+=("$TMP/$fc.bam")
    echo "   $fc: $($SAMTOOLS view -c "$TMP/$fc.bam") reads"
  done
  $SAMTOOLS merge -f -o "$TMP/merged.bam" "${parts[@]}"
  $SAMTOOLS sort -o hg001_ont_regions.bam "$TMP/merged.bam"
  $SAMTOOLS index hg001_ont_regions.bam
  fi
  echo "   hg001_ont_regions.bam: $($SAMTOOLS view -c hg001_ont_regions.bam) reads"

  echo "== longreads: trio genotypes (NA12878, father NA12891, mother NA12892), 1000 Genomes 30x"
  K=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
  vcfs=()
  for chr in chr11 chr15 chrX chr21; do
    f=1kGP_high_coverage_Illumina.$chr.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    [ $chr = chrX ] && f=1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
    regs=$(regions_bed | awk -v c=$chr '$1 == c { printf "%s%s:%d-%d", (n++ ? "," : ""), $1, $2 + 1, $3 }')
    mkdir -p "$TMP/v$chr"
    ( cd "$TMP/v$chr" && bcftools view -s NA12878,NA12891,NA12892 -r "$regs" -v snps -m2 -M2 \
        -Oz -o "$TMP/$chr.vcf.gz" "$K/$f" 2> /dev/null )
    vcfs+=("$TMP/$chr.vcf.gz")
  done
  bcftools concat "${vcfs[@]}" 2> /dev/null | bcftools annotate -x INFO -Oz -o hg001_trio_regions.vcf.gz
  tabix -f -p vcf hg001_trio_regions.vcf.gz
  echo "   $(bcftools view -H hg001_trio_regions.vcf.gz | wc -l) biallelic SNPs"

  echo "== longreads: GM12878 WGBS (ENCODE ENCSR890UQO, bed9+: rep 1 ENCFF570TIL, rep 2 ENCFF614QHA) in the same regions"
  for x in ENCFF570TIL:gm12878_wgbs_regions.bed.gz ENCFF614QHA:gm12878_wgbs_rep2_regions.bed.gz; do
    acc=${x%%:*}; out=${x##*:}
    [ -s $out ] || curl -sL https://www.encodeproject.org/files/$acc/@@download/$acc.bed.gz | gzip -dc \
      | awk 'BEGIN { OFS = "\t" } NR == FNR { c[NR] = $1; s[NR] = $2; e[NR] = $3; n = NR; next }
             { for (i = 1; i <= n; i++) if ($1 == c[i] && $2 >= s[i] && $3 <= e[i]) { print; break } }' \
           longread_regions.bed - | gzip > $out
    echo "   $out: $(gzip -dc $out | wc -l) CpG rows"
  done

  echo "== longreads: ONT haplotype-resolved 5mC pileups of HG001 (wf-human-variation), chrX and chr20"
  W=https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/wf-human-variation/sup/HG001/PAW79146/output
  for hp in 1 2; do
    f=hg001_hp${hp}_chrX_chr20_5mc.bed.gz
    [ -s $f ] && continue
    curl -s --retry 5 $W/SAMPLE.wf_mods.$hp.bedmethyl.gz | gzip -dc \
      | awk 'BEGIN { OFS = "\t" } ($1 == "chrX" || $1 == "chr20") && $4 == "m" { print $1, $2, $3, $10, $12 }' \
      | gzip > $f
    echo "   $f: $(gzip -dc $f | wc -l) CpGs (chrom, start, end, valid calls, 5mC calls)"
  done
  [ -s hg38_cpgIslandExt.txt.gz ] || curl -s -o hg38_cpgIslandExt.txt.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz

  echo "== longreads: EcoGII (m6A) chromatin-accessibility reads, HG002, and the untreated control"
  C=https://ont-open-data.s3.amazonaws.com/chrom_acc_2025.06/basecalls
  for x in PAY22766:hg002_ecogii PBA15131:hg002_native; do
    fc=${x%%:*}; name=${x##*:}
    [ -s ${name}_regions.bam.bai ] && { echo "   ${name}_regions.bam exists"; continue; }
    mkdir -p "$TMP/c$fc"
    ( cd "$TMP/c$fc" && $SAMTOOLS view -b -o "$TMP/$name.bam" "$C/$fc/basecalls.bam" \
        chr21:33000000-34000000 2> /dev/null ) || true
    $SAMTOOLS sort -o ${name}_regions.bam "$TMP/$name.bam"
    $SAMTOOLS index ${name}_regions.bam
    echo "   ${name}_regions.bam: $($SAMTOOLS view -c ${name}_regions.bam) reads"
  done
}

case $what in
  arrays) arrays ;;
  clocks) clocks ;;
  longreads) longreads ;;
  all) arrays; clocks; longreads ;;
  *) echo "usage: $0 [arrays|clocks|longreads|all]"; exit 1 ;;
esac
