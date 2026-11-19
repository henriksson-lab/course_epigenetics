"""Extract the reads of one region from a remote, coordinate-sorted BAM file
using HTTP range requests, without downloading the whole file and without an
index. Used by prepare-meth.sh.

BAM files are BGZF: a series of independently compressed blocks of at most
64 kB. So we can jump to any byte offset, find the next block, decompress it,
find the first complete alignment record and read its position. A binary
search over byte offsets then finds the part of the file that holds the
region, and only that part is downloaded.

usage: python3 -I bam_slice_remote.py URL CHROM START END OUT.bam
  START, END: 1-based, inclusive (reads whose leftmost position is in the range)
  OUT.bam is written uncompressed; pipe it through `samtools view -b` if needed.
"""
import sys, struct, zlib, urllib.request
import concurrent.futures as cf

url, chrom, start, end, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]


def fetch(a, b):
    req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
    err = None
    for _ in range(5):
        try:
            return urllib.request.urlopen(req, timeout=300).read()
        except Exception as e:  # retry on network errors
            err = e
    raise err


def total_size():
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    r = urllib.request.urlopen(req, timeout=120)
    return int(r.headers["Content-Range"].split("/")[1])


def blocks(buf):
    """Yield (offset, decompressed bytes) for each complete BGZF block in buf."""
    i = 0
    while True:
        j = buf.find(b"\x1f\x8b\x08\x04", i)
        if j < 0 or j + 18 > len(buf):
            return
        if struct.unpack_from("<H", buf, j + 10)[0] == 6 and buf[j + 12:j + 14] == b"BC":
            bsize = struct.unpack_from("<H", buf, j + 16)[0] + 1
            if j + bsize > len(buf):
                return
            try:
                yield j, zlib.decompress(buf[j + 18:j + bsize - 8], -15)
                i = j + bsize
                continue
            except zlib.error:
                pass
        i = j + 1


def valid(d, p, nref):
    """Does a plausible BAM alignment record start at byte p of d?"""
    if p + 36 > len(d):
        return False
    bs, ref, pos, lrn = struct.unpack_from("<iiiB", d, p)
    if not (32 < bs < 100000 and -1 <= ref < nref and -1 <= pos < 3_000_000_000 and lrn >= 2):
        return False
    name = d[p + 36:p + 36 + lrn]
    return len(name) == lrn and name[-1] == 0 and all(33 <= c < 127 for c in name[:-1])


def first_record(d, nref):
    """Offset of the first record in d, confirmed by the two records after it."""
    for p in range(min(len(d), 70000)):
        if valid(d, p, nref):
            q = p + 4 + struct.unpack_from("<i", d, p)[0]
            if valid(d, q, nref):
                q2 = q + 4 + struct.unpack_from("<i", d, q)[0]
                if q2 >= len(d) - 40 or valid(d, q2, nref):
                    return p
    return None


# header: magic, text, reference names
hd = b"".join(d for _, d in blocks(fetch(0, 4_000_000)))
assert hd[:4] == b"BAM\1", "not a BAM file"
p = 8 + struct.unpack_from("<i", hd, 4)[0]
nref = struct.unpack_from("<i", hd, p)[0]
p += 4
refs = []
for _ in range(nref):
    ln = struct.unpack_from("<i", hd, p)[0]
    refs.append(hd[p + 4:p + 4 + ln - 1].decode())
    p += 4 + ln + 4
header, tid = hd[:p], refs.index(chrom)
N = total_size()


def key_at(off):
    """(tid, pos) of the first record in the block after byte off; unmapped sorts last."""
    ds = list(blocks(fetch(off, off + 400_000)))
    for k in range(len(ds) - 1):
        d = ds[k][1] + ds[k + 1][1]
        r = first_record(d, nref)
        if r is not None and r < len(ds[k][1]):
            t, pos = struct.unpack_from("<ii", d, r + 4)
            return (t if t >= 0 else 10 ** 6, pos)
    return None


def search(target):
    """A byte offset before the first record >= target (within 2 MB)."""
    lo, hi = 0, N
    while hi - lo > 2_000_000:
        mid = (lo + hi) // 2
        k = key_at(mid)
        if k is not None and k < target:
            lo = mid
        else:
            hi = mid
    return lo


a = search((tid, start - 1))
b = min(search((tid, end)) + 4_000_000, N - 1)
print(f"{refs[tid]}:{start}-{end}: bytes {a}-{b} ({(b - a) / 1e6:.1f} MB of {N / 1e9:.1f} GB)", file=sys.stderr)
step = 16_000_000
parts = [(s, min(s + step - 1, b)) for s in range(a, b + 1, step)]
with cf.ThreadPoolExecutor(6) as ex:
    data = b"".join(ex.map(lambda r: fetch(*r), parts))
raw = b"".join(d for _, d in blocks(data))
p, n, recs = first_record(raw, nref), 0, bytearray()
while p + 4 <= len(raw):
    bs = struct.unpack_from("<i", raw, p)[0]
    if p + 4 + bs > len(raw):
        break
    t, pos = struct.unpack_from("<ii", raw, p + 4)
    key = (t if t >= 0 else 10 ** 6, pos)
    if key > (tid, end - 1):
        break
    if key >= (tid, start - 1):
        recs += raw[p:p + 4 + bs]
        n += 1
    p += 4 + bs
print(f"records: {n}", file=sys.stderr)
with open(out, "wb") as f:
    f.write(header)
    f.write(recs)
