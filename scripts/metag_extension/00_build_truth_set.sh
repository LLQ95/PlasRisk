#!/usr/bin/env bash
# 00_build_truth_set.sh
# ---------------------------------------------------------------------------
# Build the closed-plasmid truth set required by 01_pilot_fragment_grid.sh:
#   truth/closed_plasmids.fasta        (complete, circular reference plasmids)
#   truth/closed_plasmids_truth.tsv    (plasmid_id, length_bp, S_whole, grade,
#                                      S_full, grade_full, high_risk_genes,
#                                      cargo_genes)
#
# Source: PLSDB 2025 (data version 2024_05_31_v2), Figshare article 27252609,
# DOI 10.6084/m9.figshare.27252609.v2 (NAR 2025;53:D189-D197). PLSDB archives
# complete plasmid records; the FASTA archive is the full nucleotide set.
#
# Stage A (stdlib only): download + extract + deterministic log-size
#   stratified subsample -> closed_plasmids.fasta
# Stage B (needs plasrisk): annotate + score intact plasmids (whole sequence)
#   -> closed_plasmids_truth.tsv
#
# Usage:
#   N_TRUTH=300 bash 00_build_truth_set.sh [TRUTH_DIR] [SEED]
# Env:
#   N_TRUTH   number of reference plasmids to keep (default 300)
#   MIN_LEN   minimum plasmid length, bp (default 2000)
#   MAX_LEN   maximum plasmid length, bp (default 400000)
#   THREADS   threads for annotation (default 8)
# ---------------------------------------------------------------------------
set -euo pipefail

TRUTH_DIR="${1:-truth}"
SEED="${2:-42}"
N_TRUTH="${N_TRUTH:-300}"
MIN_LEN="${MIN_LEN:-2000}"
MAX_LEN="${MAX_LEN:-400000}"
THREADS="${THREADS:-8}"
FIGSHARE_ARTICLE="27252609"

mkdir -p "${TRUTH_DIR}" dl
CLOSED="${TRUTH_DIR}/closed_plasmids.fasta"
TRUTH_TSV="${TRUTH_DIR}/closed_plasmids_truth.tsv"

echo "=== Stage A: download PLSDB FASTA archive and subsample ==="
python3 - "${FIGSHARE_ARTICLE}" "${CLOSED}" "${SEED}" \
        "${N_TRUTH}" "${MIN_LEN}" "${MAX_LEN}" <<'PYEOF'
import glob, io, json, math, os, random, shutil, sys, tarfile, urllib.request, zipfile

article, out_fasta, seed, n_target, min_len, max_len = (
    sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]),
    int(sys.argv[5]), int(sys.argv[6]))

dl_dir = "dl"
os.makedirs(dl_dir, exist_ok=True)

# --- 1) resolve the FASTA archive download URL via the Figshare API --------
# Figshare returns HTTP 403 to the default Python-urllib User-Agent, so the
# requests identify as curl (same approach used by make_accessions.sh).
HEADERS = {"User-Agent": "curl/8.4.0", "Accept": "application/json"}


def open_url(url, timeout=120):
    req = urllib.request.Request(url, headers=HEADERS)
    return urllib.request.urlopen(req, timeout=timeout)


api = "https://api.figshare.com/v2/articles/%s/files" % article
with open_url(api, 60) as fh:
    files = json.load(fh)
fasta_files = [f for f in files
               if any(k in f["name"].lower()
                      for k in ("fasta", "fna", "nucleotide", "sequence"))
               and f["name"].lower().endswith(
                   (".tar.gz", ".tgz", ".gz", ".zip"))]
if not fasta_files:
    # fall back to the largest archive
    fasta_files = sorted(files, key=lambda f: f.get("size", 0), reverse=True)[:1]
chosen = sorted(fasta_files, key=lambda f: f.get("size", 0), reverse=True)[0]
print("  archive: %s (%.1f MB)" % (chosen["name"],
                                  chosen.get("size", 0) / 1e6))

archive = os.path.join(dl_dir, chosen["name"])
if not os.path.exists(archive):
    print("  downloading ...")
    with open_url(chosen["download_url"], 600) as r, \
            open(archive, "wb") as out:
        shutil.copyfileobj(r, out, length=1024 * 1024)
else:
    print("  cached: %s" % archive)

# --- 2) extract ------------------------------------------------------------
ext_dir = os.path.join(dl_dir, "extracted")
os.makedirs(ext_dir, exist_ok=True)
if archive.lower().endswith(".zip"):
    with zipfile.ZipFile(archive) as zf:
        names = [n for n in zf.namelist()
                 if n.lower().endswith((".fna", ".fasta", ".fa", ".fsa"))]
        if not names:  # extract everything and locate fasta afterwards
            zf.extractall(ext_dir)
        else:
            for n in names:
                zf.extract(n, ext_dir)
else:
    opener, mode = (tarfile.open, "r:gz") \
        if archive.endswith((".tar.gz", ".tgz")) else (tarfile.open, "r:*")
    with opener(archive, mode) as tf:
        members = [m for m in tf.getmembers()
                   if m.isfile() and m.name.lower().endswith(
                       (".fna", ".fasta", ".fa", ".fsa"))]
        if not members:  # extract everything and locate fasta afterwards
            tf.extractall(ext_dir)
        else:
            for m in members:
                tf.extract(m, ext_dir)

fasta_paths = glob.glob(os.path.join(ext_dir, "**", "*.f*a"), recursive=True)
fasta_paths = [p for p in fasta_paths
               if p.lower().endswith((".fna", ".fasta", ".fa", ".fsa"))]
if not fasta_paths:
    sys.exit("No FASTA found in PLSDB archive")
print("  fasta file(s): %d" % len(fasta_paths))

# --- 3) index records (path, byte offset, length) without loading all ------
records = []
for fp in fasta_paths:
    cur_id, cur_off, cur_len = None, None, 0
    with open(fp, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if cur_id is not None:
                    records.append((fp, cur_id, cur_off, cur_len))
                cur_id = line[1:].split()[0].decode()
                cur_off = fh.tell() - len(line)
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_id is not None:
            records.append((fp, cur_id, cur_off, cur_len))
print("  total PLSDB records indexed: %d" % len(records))

elig = [r for r in records if min_len <= r[3] <= max_len]
print("  within %d-%d bp: %d" % (min_len, max_len, len(elig)))
if len(elig) < n_target:
    n_target = len(elig)

# --- 4) deterministic log-size stratified subsample -----------------------
rng = random.Random(seed)
nbins = 10
edges = [min_len * (max_len / min_len) ** (i / nbins) for i in range(nbins + 1)]
bins = [[] for _ in range(nbins)]
for r in elig:
    idx = min(nbins - 1,
              max(0, int((math.log10(r[3]) - math.log10(min_len))
                         / (math.log10(max_len) - math.log10(min_len)) * nbins)))
    bins[idx].append(r)
# equal quota per non-empty size bin; fill remainder proportionally
nonempty = [b for b in bins if b]
quota = n_target // len(nonempty)
chosen_records = []
for b in nonempty:
    rng.shuffle(b)
    chosen_records.extend(b[:quota])
remaining = n_target - len(chosen_records)
pool = [r for b in nonempty for r in b[quota:]]
rng.shuffle(pool)
chosen_records.extend(pool[:max(0, remaining)])
chosen_records = chosen_records[:n_target]
chosen_ids = {r[1] for r in chosen_records}
print("  selected %d plasmids across %d size bins" %
      (len(chosen_records), len(nonempty)))

# --- 5) write closed_plasmids.fasta by byte-range copy ---------------------
n_written = 0
with open(out_fasta, "wb") as out:
    for fp in fasta_paths:
        file_recs = sorted(
            [(cid, off, length) for f, cid, off, length in records if f == fp],
            key=lambda x: x[1])
        for i, (cid, start, length) in enumerate(file_recs):
            if cid not in chosen_ids:
                continue
            end = (file_recs[i + 1][1]
                   if i + 1 < len(file_recs) else None)
            with open(fp, "rb") as src:
                src.seek(start)
                block = src.read() if end is None else src.read(end - start)
            out.write(block if block.endswith(b"\n") else block + b"\n")
            n_written += 1
print("  wrote %s (%d sequences)" % (out_fasta, n_written))
PYEOF

echo "=== Stage B: annotate + score intact plasmids (whole sequence) ==="
python3 - "${CLOSED}" "${TRUTH_TSV}" "${THREADS}" <<'PYEOF'
import sys
from plasrisk import annotate_fasta, PlasRiskScorer

fasta, out_tsv, threads = sys.argv[1], sys.argv[2], int(sys.argv[3])

ann = annotate_fasta(fasta, use_abricate=True, threads=threads)
feats = {f.seq_id: f for f in ann.features}
if ann.warnings:
    for w in ann.warnings:
        print("  warning:", w)

full = PlasRiskScorer(mode="full")
lite = PlasRiskScorer(mode="lite")

cols = ["plasmid_id", "length_bp", "S_whole", "grade", "S_full",
        "grade_full", "high_risk_genes", "cargo_genes"]
rows = []
for f in ann.features:
    rf = full.score(f)
    rl = lite.score(f)
    cargo = []
    for name in list(f.arg_names) + list(f.bm_gene_names) + list(f.vf_names):
        if name not in cargo:
            cargo.append(name)
    rows.append({
        "plasmid_id": f.seq_id,
        "length_bp": f.length_bp,
        "S_whole": rl["S_norm"],          # target a perfect bin converges to
        "grade": rl["grade"],
        "S_full": rf["S_norm"],
        "grade_full": rf["grade"],
        "high_risk_genes": rf.get("high_risk_genes", ""),
        "cargo_genes": ";".join(cargo),
    })

rows.sort(key=lambda r: r["plasmid_id"])
with open(out_tsv, "w", encoding="utf-8", newline="") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print("  wrote %s (%d plasmids)" % (out_tsv, len(rows)))
PYEOF

echo "DONE. Truth set ready in ${TRUTH_DIR}/"
