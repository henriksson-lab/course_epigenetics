#!/bin/bash
# Download the course datasets listed in data/manifest.tsv into data/.
# The files are large and never checked in (see .gitignore).
#
#   data/download.sh                    # everything
#   data/download.sh genomic-intervals  # only files used by one unit
#
# Already downloaded files are skipped. Works the same on a laptop or in a
# SLURM job (sbatch --wrap "data/download.sh").
set -euo pipefail
cd "$(dirname "$0")"
unit="${1:-}"
grep -v '^#' manifest.tsv | while IFS=$'\t' read -r id url desc used_by; do
  [ -n "$unit" ] && [[ "$used_by" != *"$unit"* ]] && continue
  if [ -s "$id" ]; then echo "have $id"; continue; fi
  echo "get  $id  ($desc)"
  curl -fL -sS --retry 3 -o "$id.part" "$url" && mv "$id.part" "$id"
done
