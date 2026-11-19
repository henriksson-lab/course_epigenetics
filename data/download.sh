#!/bin/bash
# Get the course datasets listed in data/manifest.tsv into data/.
# The files are large and never checked in (see .gitignore).
#
#   data/download.sh                    # everything
#   data/download.sh genomic-intervals  # only files used by one unit
#
# The url column of the manifest is one of:
#   https://...        downloaded with curl
#   prepare-*.sh       made by that script in data/ (e.g. prepare-reads.sh: ENCODE
#                      reads subset to chr21+chrM); each script is run once
#   derived            built here from other files (bowtie2 index, bgzipped genome)
# Files that already exist are skipped. Works the same on a laptop or in a
# SLURM job (sbatch --wrap "data/download.sh <unit>").
set -euo pipefail
cd "$(dirname "$0")"
unit="${1:-}"
prepare=()
derived=()
while IFS=$'\t' read -r id url desc used_by; do
  [[ -z "$id" || "$id" == \#* ]] && continue
  [ -n "$unit" ] && [[ "$used_by" != *"$unit"* ]] && continue
  if [ -s "$id" ] || [ -s "$id.1.bt2" ]; then echo "have $id"; continue; fi
  case "$url" in
    http*)            echo "get  $id  ($desc)"
                      curl -fL -sS --retry 3 -o "$id.part" "$url" && mv "$id.part" "$id" ;;
    prepare-*.sh)     [[ " ${prepare[*]-} " == *" $url "* ]] || prepare+=("$url") ;;
    derived)          derived+=("$id") ;;
  esac
done < manifest.tsv

for script in "${prepare[@]+"${prepare[@]}"}"; do
  echo "run $script (builds files that are not plain downloads; see its header)"
  ./"$script"
done

for id in "${derived[@]+"${derived[@]}"}"; do
  case "$id" in
    bowtie2/hg38_chr21_chrM)
      [ -s hg38_chr21_chrM.fa ] || ./prepare-reads.sh
      echo "build bowtie2 index $id"
      mkdir -p bowtie2 && bowtie2-build --threads 4 hg38_chr21_chrM.fa bowtie2/hg38_chr21_chrM > /dev/null ;;
    hg38.fa.bgz)
      [ -s hg38.fa.gz ] || curl -fL -sS --retry 3 -o hg38.fa.gz \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
      echo "bgzip and index $id"
      zcat hg38.fa.gz | bgzip -@ 4 > hg38.fa.bgz && samtools faidx hg38.fa.bgz ;;
  esac
done
