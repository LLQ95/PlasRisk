#!/usr/bin/env bash
# 00_build_truth_set.sh
# ---------------------------------------------------------------------------
# Build the closed-plasmid truth set required by 01_pilot_fragment_grid.sh:
#   truth/closed_plasmids.fasta        (complete reference plasmids)
#   truth/closed_plasmids_truth.tsv    (plasmid_id, length_bp, S_whole, grade,
#                                      S_full, grade_full, high_risk_genes,
#                                      cargo_genes)
#
# Source: PLSDB 2025 (NAR 2025;53:D189-D197). The complete nucleotide set is
# streamed from the official PLSDB server (Saarland, DE):
#     https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta
# This endpoint returns an uncompressed sequences.fasta (~7.4 GB) and does NOT
# go through Figshare. The Figshare API (api.figshare.com) is unreachable from
# some networks (HTTP 403 at the AWS ELB), so it is only used as a fallback.
#
# Stage A (stdlib only): stream the FASTA once, reservoir-sample within
#   prespecified log-size bins -> closed_plasmids.fasta (the full 7.4 GB is
#   read but discarded; it is never written to disk).
# Stage B (needs plasrisk): annotate + score intact plasmids (whole sequence)
#   -> closed_plasmids_truth.tsv
#
# Usage:
#   N_TRUTH=300 bash 00_build_truth_set.sh [TRUTH_DIR] [SEED]
# Env:
#   N_TRUTH     number of reference plasmids to keep (default 300)
#   MIN_LEN     minimum plasmid length, bp (default 2000)
#   MAX_LEN     maximum plasmid length, bp (default 400000)
#   THREADS     threads for annotation (default 8)
#   PLSDB_FASTA  use a pre-downloaded sequences.fasta instead of streaming
#               (download once, resumable:
#                curl -L -C - -A "curl/8.4.0" \
#                  https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta \
#                  -o dl/sequences.fasta)
# ---------------------------------------------------------------------------
set -euo pipefail

TRUTH_DIR="${1:-truth}"
SEED="${2:-42}"
N_TRUTH="${N_TRUTH:-300}"
MIN_LEN="${MIN_LEN:-2000}"
MAX_LEN="${MAX_LEN:-400000}"
THREADS="${THREADS:-8}"
PLSDB_FASTA="${PLSDB_FASTA:-}"
PLSDB_DIRECT="https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta"

mkdir -p "${TRUTH_DIR}" dl
CLOSED="${TRUTH_DIR}/closed_plasmids.fasta"
TRUTH_TSV="${TRUTH_DIR}/closed_plasmids_truth.tsv"

echo "=== Stage A: stream PLSDB FASTA and reservoir-sample ==="
python3 - "${CLOSED}" "${SEED}" "${N_TRUTH}" "${MIN_LEN}" "${MAX_LEN}" \
        "${PLSDB_FASTA}" "${PLSDB_DIRECT}" <<'PYEOF'
import io, math, os, random, sys, urllib.request

(out_fasta, seed, n_target, min_len, max_len,
 local_fasta, direct_url) = (
    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]),
    int(sys.argv[5]), sys.argv[6], sys.argv[7])

nbins = 10
# One reservoir per log-size bin. Capacity >= n_target so the combined stored
# candidates always suffice even if all eligible records fall in a single bin.
cap = max(200, n_target)
HEADERS = {"User-Agent": "curl/8.4.0", "Accept": "*/*"}


def bin_index(length):
    if length < min_len or length > max_len:
        return -1
    lo, hi = math.log10(min_len), math.log10(max_len)
    idx = int((math.log10(length) - lo) / (hi - lo) * nbins)
    return min(nbins - 1, max(0, idx))


def open_source():
    """Yield raw lines from a local file or the remote PLSDB endpoint."""
    if local_fasta:
        if not os.path.exists(local_fasta):
            sys.exit("PLSDB_FASTA set but not found: %s" % local_fasta)
        print("  source: local %s" % local_fasta)
        return open(local_fasta, "rb")
    print("  source: %s" % direct_url)
    print("  (streaming ~7.4 GB; records are read and discarded, not stored)")
    req = urllib.request.Request(direct_url, headers=HEADERS)
    try:
        return urllib.request.urlopen(req, timeout=120)
    except Exception as exc:
        sys.exit(
            "Could not reach the PLSDB direct endpoint (%s).\n"
            "Pre-download with:\n"
            "  curl -L -C - -A 'curl/8.4.0' %s -o dl/sequences.fasta\n"
            "then re-run: PLSDB_FASTA=dl/sequences.fasta bash 00_build_truth_set.sh\n"
            "Underlying error: %s" % (direct_url, direct_url, exc))


# --- stream once and fill per-bin reservoirs (Algorithm R) ----------------
reservoirs = [[] for _ in range(nbins)]   # list of (header_bytes, seq_bytes)
seen = [0] * nbins
rng = random.Random(seed)

cur_header, cur_seq, cur_too_big = None, bytearray(), False
n_total = 0


def finalize(header, seq, too_big):
    global n_total
    if header is None:
        return
    n_total += 1
    if too_big:
        return
    length = len(seq)
    bi = bin_index(length)
    if bi < 0:
        return
    k = seen[bi]
    seen[bi] = k + 1
    rec = (header, bytes(seq))
    if len(reservoirs[bi]) < cap:
        reservoirs[bi].append(rec)
    else:
        j = rng.randint(0, k)
        if j < cap:
            reservoirs[bi][j] = rec


src = open_source()
try:
    for line in src:
        if line.startswith(b">"):
            finalize(cur_header, cur_seq, cur_too_big)
            cur_header, cur_seq, cur_too_big = line, bytearray(), False
        else:
            if not cur_too_big:
                cur_seq.extend(b"".join(line.split()))
                if len(cur_seq) > max_len:
                    cur_too_big = True
                    cur_seq = bytearray()
    finalize(cur_header, cur_seq, cur_too_big)
finally:
    src.close()

print("  total PLSDB records scanned: %d" % n_total)
nonempty = [i for i in range(nbins) if seen[i] > 0]
print("  eligible records by size bin: %s" %
      ", ".join("%d" % seen[i] for i in nonempty))
total_elig = sum(seen[i] for i in nonempty)
if total_elig < n_target:
    print("  WARNING: only %d eligible records; using all" % total_elig)
    n_target = total_elig

# --- deterministic stratified draw from the reservoirs --------------------
quota = n_target // len(nonempty)
chosen, used = [], set()
for bi in nonempty:
    pool = reservoirs[bi]
    take = min(quota, len(pool))
    for rec in rng.sample(pool, take):
        chosen.append(rec)
# fill remainder from all unused reservoir candidates
remaining = n_target - len(chosen)
if remaining > 0:
    rest = []
    for bi in nonempty:
        for rec in reservoirs[bi]:
            if rec not in chosen:
                rest.append(rec)
    rng.shuffle(rest)
    chosen.extend(rest[:remaining])
chosen = chosen[:n_target]

with open(out_fasta, "wb") as out:
    for header, seq in chosen:
        out.write(header if header.endswith(b"\n") else header + b"\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + b"\n")
print("  wrote %s (%d plasmids)" % (out_fasta, len(chosen)))
PYEOF

echo "=== Stage B: annotate + score intact plasmids (whole sequence) ==="
python3 - "${CLOSED}" "${TRUTH_TSV}" "${THREADS}" <<'PYEOF'
import sys
from plasrisk import annotate_fasta, PlasRiskScorer

fasta, out_tsv, threads = sys.argv[1], sys.argv[2], int(sys.argv[3])

ann = annotate_fasta(fasta, use_abricate=True, threads=threads)
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
