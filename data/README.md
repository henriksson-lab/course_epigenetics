# Course data

Nothing large is checked in. Tracked here are only:

- `manifest.tsv`: every downloadable file (what, from where, which unit uses it);
- `download.sh`: fetches them, `data/download.sh` or `data/download.sh <unit-slug>`;
- `prepare-reads.sh` + `filter_fastq_by_names.py`: build the small read sets
  for the alignment, QC and ATAC-seq units (see below).

All genomic data are **GRCh38** (the Corces 2016 counts are hg19 and are lifted
over in the exercise). ENCODE accessions are in the URLs.

## Prepared read sets (not downloadable as such)

`prepare-reads.sh` makes, from ENCODE GM12878 data:

| file | source | reads |
|---|---|---|
| `gm12878_atac_rep{1,2}_R{1,2}.fastq.gz` | ATAC-seq ENCSR095QNB, PE 2×101 | pairs ENCODE aligned to chr21 or chrM |
| `gm12878_ctcf_chip.fastq.gz` | CTCF ChIP-seq ENCSR000DZN rep 1, SE 36 | reads aligned to chr21 or chrM |
| `gm12878_ctcf_input.fastq.gz` | input ENCSR000EYX, SE 28 | reads aligned to chr21 or chrM |
| `hg38_chr21_chrM.fa` | UCSC hg38 | the matching small reference |

The reads are the original, untrimmed FASTQ records, so trimming, alignment
and filtering are real; only the choice of reads is a shortcut. Running the
script needs about 25 GB of temporary space and an hour.

**For students:** either run `data/prepare-reads.sh` once on the cluster, or
the teacher places the prepared files (a few hundred MB) in a shared folder.
**[?]** Decide where to host them.

The bowtie2 index (`data/bowtie2/`) and `hg38.fa.bgz` (bgzipped + `samtools faidx`)
are built from these files with the commands in the units.
