#!/usr/bin/env python3
"""Merge plan/manifest-<block>.tsv fragments into data/manifest.tsv and
recompute every row's used_by from what the units actually reference.

    python3 plan/merge-manifests.py

A unit "uses" a file if its lec.qmd, exercise.qmd, quiz.qmd or make-*-figs.*
mentions the file name (for derived directories such as bowtie2/, the base name).
Declared used_by entries are kept as well. Rows keep their order; new rows
from fragments are appended in fragment order, skipping ids already present.
"""
import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = os.path.join(ROOT, "data", "manifest.tsv")


def read_rows(path):
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            parts += [""] * (4 - len(parts))
        rows.append(parts[:4])
    return rows


header = [l.rstrip("\n") for l in open(MAN) if l.startswith("#")]
rows = read_rows(MAN)
by_id = {r[0]: r for r in rows}
for frag in sorted(glob.glob(os.path.join(ROOT, "plan", "manifest-*.tsv"))):
    for r in read_rows(frag):
        if r[0] in by_id:
            old = by_id[r[0]]
            units = {u.strip() for u in (old[3] + "," + r[3]).split(",") if u.strip()}
            old[3] = ", ".join(sorted(units))
        else:
            rows.append(r)
            by_id[r[0]] = r

units = sorted(d for d in os.listdir(ROOT)
               if os.path.isfile(os.path.join(ROOT, d, "lec.qmd")))
text = {}
for u in units:
    files = [os.path.join(ROOT, u, f) for f in ("lec.qmd", "exercise.qmd", "quiz.qmd")]
    files += glob.glob(os.path.join(ROOT, u, "make-*"))
    text[u] = "\n".join(open(f, errors="replace").read() for f in files if os.path.isfile(f))

for r in rows:
    name = os.path.basename(r[0].rstrip("/"))
    pat = re.compile(r"(?<![\w.-])" + re.escape(name) + r"(?![\w-])")
    found = {u for u in units if pat.search(text[u])}
    declared = {x.strip() for x in r[3].split(",") if x.strip()}
    r[3] = ", ".join(sorted(declared | found))

with open(MAN, "w") as f:
    f.write("\n".join(header) + "\n" if header else "")
    for r in rows:
        f.write("\t".join(r) + "\n")
print(f"{len(rows)} rows; units scanned: {len(units)}")
for r in rows:
    if not r[3]:
        print("  unused:", r[0])
