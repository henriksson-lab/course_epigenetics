#!/bin/bash
# Build the small TCGA methylation tables used by the unit
# "Epigenetics in disease and therapy" (cancer-epigenetics).
#
#   data/prepare-dis.sh            # about 10 minutes, ~1 GB temporary space, 4 GB RAM
#
# Sources (all public, no login):
#   - UCSC Xena TCGA hub: Illumina HumanMethylation450 beta values (TCGA legacy
#     level 3) for colon adenocarcinoma (COAD) and lower-grade glioma (LGG),
#     and the COAD clinical matrix (tumour location).
#   - Zhou, Laird & Shen 2017 (InfiniumAnnotation): 450k probe annotation on
#     GRCh38 with GENCODE v36 genes, CpG-island context and a recommended mask.
#   - cBioPortal REST API: TCGA PanCancer Atlas MSI scores (MANTIS), molecular
#     subtype, BRAF/KRAS mutations and MLH1 RNA (COAD); IDH status and
#     pan-glioma methylation cluster (Ceccarelli et al. 2016) for gliomas.
#
# Writes into data/:
#   tcga_coad_450k_beta.tsv.gz    beta values, ~41 000 probes x 337 samples
#   tcga_coad_samples.tsv         one row per COAD sample (tumour/normal, MSI, BRAF, ...)
#   tcga_lgg_450k_beta.tsv.gz     beta values, ~20 000 probes x 516 gliomas
#   tcga_lgg_samples.tsv          one row per glioma (IDH status, 1p/19q, cluster)
#   tcga_450k_probes.tsv.gz       annotation of the probes in both tables (GRCh38)
#
# Probe choice (a shortcut, said so in the unit): after removing masked,
# non-autosomal and mostly-missing probes, a random 10% of probes (seed 2026)
# plus every probe in the CpG islands of a panel of genes discussed in the
# lecture. Betas are rounded to 3 decimals.
set -euo pipefail
cd "$(dirname "$0")"
DATA=$(pwd)
TMP=$(mktemp -d -p "$DATA" tmp-dis.XXXX)
trap 'rm -rf "$TMP"' EXIT

XENA=https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download
ZHOU=https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450
CBIO=https://www.cbioportal.org/api
CRC=coadread_tcga_pan_can_atlas_2018
GLI=lgggbm_tcga_pub
get() { curl -fL -sS --retry 3 -o "$TMP/$1" "$2"; }

echo "download (about 800 MB)"
get coad_450k.gz    "$XENA/TCGA.COAD.sampleMap%2FHumanMethylation450.gz"
get lgg_450k.gz     "$XENA/TCGA.LGG.sampleMap%2FHumanMethylation450.gz"
get coad_clin.tsv   "$XENA/TCGA.COAD.sampleMap%2FCOAD_clinicalMatrix"
get anno.tsv.gz     "$ZHOU/HM450.hg38.manifest.gencode.v36.tsv.gz"
get mask.tsv.gz     "$ZHOU/HM450.hg38.mask.tsv.gz"
for a in MSI_SCORE_MANTIS MSI_SENSOR_SCORE; do
  get "crc_$a.json" "$CBIO/studies/$CRC/clinical-data?clinicalDataType=SAMPLE&attributeId=$a&pageSize=10000"
done
for a in SUBTYPE AGE SEX; do
  get "crc_$a.json" "$CBIO/studies/$CRC/clinical-data?clinicalDataType=PATIENT&attributeId=$a&pageSize=10000"
done
get crc_sequenced.json "$CBIO/sample-lists/${CRC}_sequenced"
curl -fL -sS --retry 3 -o "$TMP/crc_mut.json" -X POST -H 'Content-Type: application/json' \
  "$CBIO/molecular-profiles/${CRC}_mutations/mutations/fetch?projection=SUMMARY" \
  -d "{\"sampleListId\":\"${CRC}_all\",\"entrezGeneIds\":[673,3845]}"
curl -fL -sS --retry 3 -o "$TMP/crc_rna.json" -X POST -H 'Content-Type: application/json' \
  "$CBIO/molecular-profiles/${CRC}_rna_seq_v2_mrna/molecular-data/fetch" \
  -d "{\"sampleListId\":\"${CRC}_all\",\"entrezGeneIds\":[4292]}"
for a in IDH_STATUS IDH_CODEL_SUBTYPE PAN_GLIOMA_DNA_METHYLATION_CLUSTER; do
  get "gli_$a.json" "$CBIO/studies/$GLI/clinical-data?clinicalDataType=SAMPLE&attributeId=$a&pageSize=10000"
done
get gli_patient.json "$CBIO/studies/$GLI/clinical-data?clinicalDataType=PATIENT&attributeId=HISTOLOGICAL_DIAGNOSIS&pageSize=10000"

echo "subset and write tables (R)"
Rscript - "$TMP" "$DATA" <<'EOF'
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
args <- commandArgs(TRUE); tmp <- args[1]; out <- args[2]
f <- function(x) file.path(tmp, x)

an <- fread(f("anno.tsv.gz")); mk <- fread(f("mask.tsv.gz"))
an$masked <- an$probeID %in% mk$Probe_ID[mk$M_general %in% TRUE]
an$cgi <- fifelse(is.na(an$CGIposition), "open sea",
          fifelse(an$CGIposition == "Island", "island",
          fifelse(grepl("Shore", an$CGIposition), "shore", "shelf")))
# distance to the nearest TSS among all annotated transcripts
an$tss_dist <- vapply(strsplit(as.character(an$distToTSS), ";"), function(v) {
  v <- suppressWarnings(as.numeric(v)); if (all(is.na(v))) NA_real_ else v[which.min(abs(v))] }, 0)
panel <- c("MLH1", "CDKN2A", "MGMT", "SEPTIN9", "SFRP2", "GATA4", "RUNX3", "NEUROG1",
           "VIM", "CACNA1G", "SOCS1", "IGF2", "CRABP1", "TIMP3", "APC", "BRCA1", "GAPDH", "ACTB")
in_panel <- vapply(strsplit(an$genesUniq, ";"), function(v) any(v %in% panel), TRUE) &
            an$cgi == "island" & abs(an$tss_dist) <= 1500

read_beta <- function(file) {
  m <- fread(file); B <- as.matrix(m[, -1]); rownames(B) <- m[[1]]; B }
keep_probes <- function(B) {
  a <- an[match(rownames(B), an$probeID)]
  ok <- !a$masked & a$CpG_chrm %in% paste0("chr", 1:22) & rowMeans(is.na(B)) < 0.1
  rownames(B)[ok] }

# ---- COAD -----------------------------------------------------------------
B <- read_beta(f("coad_450k.gz"))
ok <- keep_probes(B)
set.seed(2026)
pick <- union(sample(ok, round(0.1 * length(ok))), intersect(ok, an$probeID[in_panel]))
pick <- an$probeID[an$probeID %in% pick]                           # genome order
Bc <- round(B[pick, ], 3)
st <- substr(colnames(Bc), 14, 15)
keep_s <- st %in% c("01", "11")                                     # primary tumours and normals
Bc <- Bc[, keep_s]
cat("COAD:", length(ok), "good probes;", nrow(Bc), "kept;", ncol(Bc), "samples\n")
fwrite(data.table(probe = rownames(Bc), Bc), file.path(out, "tcga_coad_450k_beta.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)

cl <- fread(f("coad_clin.tsv"))
js <- function(x) as.data.table(fromJSON(f(x)))
msi <- js("crc_MSI_SCORE_MANTIS.json"); msis <- js("crc_MSI_SENSOR_SCORE.json")
sub <- js("crc_SUBTYPE.json"); age <- js("crc_AGE.json"); sex <- js("crc_SEX.json")
sq <- fromJSON(f("crc_sequenced.json"))$sampleIds
mut <- js("crc_mut.json"); rna <- js("crc_rna.json")
s <- colnames(Bc); p <- substr(s, 1, 12)
right <- c("Cecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon")
left  <- c("Splenic Flexure", "Descending Colon", "Sigmoid Colon", "Rectosigmoid Junction")
loc <- cl$anatomic_neoplasm_subdivision[match(s, cl$sampleID)]
smp <- data.table(
  sample = s, patient = p,
  type = ifelse(substr(s, 14, 15) == "11", "normal", "tumour"),
  side = ifelse(loc %in% right, "right", ifelse(loc %in% left, "left", NA)),
  subtype = sub$value[match(p, sub$patientId)],
  msi_mantis = as.numeric(msi$value[match(s, msi$sampleId)]),
  msisensor = as.numeric(msis$value[match(s, msis$sampleId)]),
  braf = ifelse(!s %in% sq, NA, ifelse(s %in% mut[entrezGeneId == 673 & proteinChange == "V600E"]$sampleId, "V600E", "wt")),
  kras = ifelse(!s %in% sq, NA, ifelse(s %in% mut[entrezGeneId == 3845 & mutationType == "Missense_Mutation"]$sampleId, "mut", "wt")),
  mlh1_rna = round(rna$value[match(s, rna$sampleId)], 1),
  age = as.numeric(age$value[match(p, age$patientId)]),
  sex = sex$value[match(p, sex$patientId)])
smp$msi <- ifelse(smp$type == "normal" | is.na(smp$msi_mantis), NA, ifelse(smp$msi_mantis > 0.4, "MSI-H", "MSS"))
smp[type == "normal", c("subtype", "braf", "kras", "msi_mantis", "msisensor") := NA]
fwrite(smp, file.path(out, "tcga_coad_samples.tsv"), sep = "\t", na = "NA", quote = FALSE)

# ---- LGG ------------------------------------------------------------------
G <- read_beta(f("lgg_450k.gz"))
G <- G[, substr(colnames(G), 14, 15) == "01"]
idh <- js("gli_IDH_STATUS.json"); cod <- js("gli_IDH_CODEL_SUBTYPE.json")
clu <- js("gli_PAN_GLIOMA_DNA_METHYLATION_CLUSTER.json"); his <- js("gli_patient.json")
G <- G[, colnames(G) %in% idh$sampleId]
okg <- keep_probes(G)
set.seed(2026)
pg <- union(sample(okg, 20000), intersect(okg, an$probeID[in_panel]))
pg <- an$probeID[an$probeID %in% pg]
Gs <- round(G[pg, ], 3)
cat("LGG:", length(okg), "good probes;", nrow(Gs), "kept;", ncol(Gs), "gliomas\n")
fwrite(data.table(probe = rownames(Gs), Gs), file.path(out, "tcga_lgg_450k_beta.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)
gs <- colnames(Gs)
fwrite(data.table(sample = gs,
                  idh = idh$value[match(gs, idh$sampleId)],
                  idh_codel = cod$value[match(gs, cod$sampleId)],
                  meth_cluster = clu$value[match(gs, clu$sampleId)],
                  histology = his$value[match(substr(gs, 1, 12), his$patientId)]),
       file.path(out, "tcga_lgg_samples.tsv"), sep = "\t", na = "NA", quote = FALSE)

# ---- probe annotation -------------------------------------------------------
ids <- union(rownames(Bc), rownames(Gs))
a <- an[match(ids, an$probeID)]
a <- a[order(match(CpG_chrm, paste0("chr", 1:22)), CpG_beg)]
fwrite(a[, .(probe = probeID, chr = CpG_chrm, pos = CpG_beg + 1L, genes = genesUniq,
             tss_dist, cgi)], file.path(out, "tcga_450k_probes.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)
cat("probes annotated:", nrow(a), "\n")
EOF
ls -l "$DATA"/tcga_*
