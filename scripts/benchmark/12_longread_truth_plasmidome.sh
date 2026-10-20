#!/usr/bin/env bash
# 12_longread_truth_plasmidome.sh
#
# P0-1 (plan A): technology-independent truth set from long-read
# plasmidome assemblies.
#
# Source: Kirstahler et al. 2021, mSystems 6:e00283-21,
# "A Peek into the Plasmidome of Global Sewage" (PRJEB41171), 24
# plasmid-enriched (phi29 MDA) sewage samples, 22 countries, 5
# continents. Each sample has (i) a polished single-Nanopore-read
# assembly (miniasm + racon + 5 rounds pilon with the matched Illumina
# reads), deposited as ERZ1694234-ERZ1694257, containing 159,322
# circular molecules out of 165,302 contigs, and (ii) matched Illumina
# NextSeq 2x150 reads (ERX5299122-ERX5299145).
#
# Design:
#   truth side  = circular molecules in the 24 long-read assemblies
#   test side   = matched Illumina reads re-assembled with MEGAHIT
#                 (the frozen short-read assembly path), run through
#                 geNomad end-to-end + plasrisk meta with the exact
#                 settings used in the prospective cohorts
#   comparison  = minimap2 recovery of truth molecules by the
#                 short-read assemblies, stratified by length and gene
#                 content; recovery funnel (assembled -> geNomad
#                 called -> PlasRisk graded); direct PlasRisk grading of
#                 the long-read truth molecules themselves
#
# Usage:
#   bash 12_longread_truth_plasmidome.sh \
#        /path/to/assemblies_plasmidome.tsv \
#        results/plasmidome_truth 16 2
#
# Resume-safe: every per-sample stage is skipped when its marker output
# exists. Delete the relevant output directory to rerun a stage.
#
# Outputs:
#   outdir/00_manifest/short_read_manifest.tsv
#   outdir/01_truth_longread/<id>.fna               (ERZ assemblies)
#   outdir/02_short_reads/<id>/{R1,R2}.fastq.gz
#   outdir/03_megahit/<id>/final.contigs.fa
#   outdir/04_genomad_short/<id>/, 04_genomad_long/<id>/
#   outdir/05_plasrisk_short/<id>/, 05_plasrisk_long/<id>/
#   outdir/06_circularity/<id>.circular.tsv
#   outdir/07_align/<id>.paf
#   outdir/truth_funnel.tsv         (one row per sample)
#   outdir/truth_recovery_contigs.tsv (one row per truth molecule)
#   outdir/truth_long_grades.tsv    (PlasRisk grades of truth molecules)
#   outdir/truth_summary.txt

set -uo pipefail

# English locale for all log output (timestamps, tool messages).
export LC_ALL=C

TRUTH_MANIFEST="${1:-}"
OUTDIR="${2:-results/plasmidome_truth}"
THREADS="${3:-16}"
JOBS="${4:-2}"
# minimap2 recovery thresholds
OVL_IDY="${OVL_IDY:-0.90}"
OVL_COV_STRICT="${OVL_COV_STRICT:-0.80}"
OVL_COV_RELAXED="${OVL_COV_RELAXED:-0.50}"
CIRC_MINBP="${CIRC_MINBP:-200}"
MIN_CONTIG="${MIN_CONTIG:-1000}"

if [[ -z "$TRUTH_MANIFEST" || ! -s "$TRUTH_MANIFEST" ]]; then
  echo "usage: bash 12_longread_truth_plasmidome.sh <assemblies_plasmidome.tsv> <outdir> [threads=16] [jobs=2]"
  exit 1
fi
for bin in python3 aria2c megahit minimap2 genomad plasrisk; do
  command -v "$bin" >/dev/null 2>&1 || { echo "[error] $bin not in PATH"; exit 1; }
done
GENOMAD_DB="${GENOMAD_DB:-}"
if [[ -z "$GENOMAD_DB" ]]; then
  echo "[error] GENOMAD_DB not set; export it to the geNomad database directory"
  exit 1
fi

BASEDIR="$(cd "$(dirname "$OUTDIR")" 2>/dev/null && pwd)/$(basename "$OUTDIR")"
mkdir -p "$BASEDIR"/{00_manifest,01_truth_longread,02_short_reads,03_megahit,04_genomad_short,04_genomad_long,05_plasrisk_short,05_plasrisk_long,06_circularity,07_align,logs}
LOG_DIR="$BASEDIR/logs"
FAIL_LOG="$BASEDIR/failed_samples.tsv"
: > "${FAIL_LOG}.tmp"

GE_VER="$(genomad --version 2>/dev/null | grep -oE '^[0-9]+\.[0-9]+' | head -1 || true)"
GE_ARG=""
case "$GE_VER" in
  1.12*|1.13*|1.14*|1.15*|1.16*|1.17*|1.18*|1.19*|1.2*|1.[2-9][0-9]*) GE_ARG="--genomad-version v1.12" ;;
esac
echo "[truth] threads/job: ${THREADS} | jobs: ${JOBS} | genomad: ${GE_VER:-?} | idy ${OVL_IDY} | cov ${OVL_COV_STRICT}/${OVL_COV_RELAXED}"

# ---------------- stage 0: ENA run manifest for matched Illumina reads ---------
SHORT_MANIFEST="$BASEDIR/00_manifest/short_read_manifest.tsv"
if [[ ! -s "$SHORT_MANIFEST" ]]; then
  echo "[truth] querying ENA for Illumina runs under PRJEB41171 ..."
  curl -sS --retry 3 --retry-delay 5 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB41171&result=read_run&fields=run_accession,experiment_accession,sample_alias,instrument_model,fastq_ftp,fastq_bytes&format=tsv&limit=0" \
    -o "$BASEDIR/00_manifest/ena_runs.tsv"
  python3 - "$BASEDIR/00_manifest/ena_runs.tsv" "$TRUTH_MANIFEST" "$SHORT_MANIFEST" <<'PYEOF'
import sys, os
ena, truth, out = sys.argv[1:4]
# long-read manifest: 4 cols no header: id, habitat, pair, url
id2alias = {}
with open(truth, encoding="utf-8") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[0].startswith("GS_"):
            id2alias[f[0]] = f[0]
rows = []
with open(ena, encoding="utf-8") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(header):
            continue
        instr = f[idx["instrument_model"]].lower()
        ftp = f[idx["fastq_ftp"]]
        files = [x for x in ftp.split(";") if x.strip()]
        if "nextseq" not in instr or len(files) != 2:
            continue
        alias = f[idx["sample_alias"]].strip()
        # aliases are like GS_Plasmidome_35; normalize stray spaces
        alias = alias.replace(" ", "")
        rows.append((alias, f[idx["run_accession"]],
                     f[idx["experiment_accession"]],
                     "https://" + files[0], "https://" + files[1]))
rows.sort()
missing = sorted(set(id2alias) - {r[0] for r in rows})
if missing:
    print("[warn] %d long-read assemblies have no matched Illumina run: %s"
          % (len(missing), ",".join(missing)), file=sys.stderr)
with open(out, "w", encoding="utf-8") as fh:
    fh.write("id\trun_accession\texperiment_accession\tr1_url\tr2_url\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
print("[truth] short-read manifest: %d samples" % len(rows))
PYEOF
fi
[[ -s "$SHORT_MANIFEST" ]] || { echo "[error] failed to build short-read manifest"; exit 1; }

# optional: author's per-contig circularity tables (figshare Table C)
AUTHOR_TABLES="$BASEDIR/00_manifest/author_tables"
mkdir -p "$AUTHOR_TABLES"
if [[ ! -s "$AUTHOR_TABLES/.done" ]]; then
  python3 - "$AUTHOR_TABLES" <<'PYEOF' || true
import json, os, sys, urllib.request
dest = sys.argv[1]
try:
    req = urllib.request.Request(
        "https://api.figshare.com/v2/projects/94448/articles",
        headers={"User-Agent": "plasrisk-truth"})
    arts = json.load(urllib.request.urlopen(req, timeout=60))
    for a in arts:
        title = a.get("title", "")
        if "Table" not in title and "table" not in title:
            continue
        for f in a.get("files", []):
            name = f.get("name", "")
            url = f.get("download_url", "")
            if name.startswith("Table_C") or "Table_C" in name:
                out = os.path.join(dest, name)
                urllib.request.urlretrieve(url, out)
                print("downloaded", name)
    open(os.path.join(dest, ".done"), "w").write("ok\n")
except Exception as exc:  # network may be blocked; end-overlap fallback is used
    print("[warn] figshare author tables unavailable (%s); "
          "using end-overlap circularity" % exc)
PYEOF
fi

# ---------------- per-sample stages --------------------------------------------
download_url() { # gz url dir
  local gz="$1" url="$2" d="$3"
  if [[ -s "$gz" ]] && gzip -t "$gz" 2>/dev/null; then return 0; fi
  rm -f "$gz"
  aria2c -x 4 -s 4 -k 1M --file-allocation=none --console-log-level=warn \
    --summary-interval=0 --retry-wait=5 -m 3 -d "$d" \
    -o "$(basename "$gz")" "$url" || return 1
  gzip -t "$gz" 2>/dev/null || return 1
}

process_one() {
  local sid="$1" url_long="$2" r1="$3" r2="$4"
  local lg="$LOG_DIR/${sid}.log"
  local long_gz="$BASEDIR/01_truth_longread/${sid}.fna.gz"
  local long_fa="$BASEDIR/01_truth_longread/${sid}.fna"
  local rdir="$BASEDIR/02_short_reads/$sid"
  local r1gz="$rdir/R1.fastq.gz" r2gz="$rdir/R2.fastq.gz"
  local mh="$BASEDIR/03_megahit/$sid"
  local short_fa="$mh/final.contigs.fa"
  local gs="$BASEDIR/04_genomad_short/$sid"
  local gl="$BASEDIR/04_genomad_long/$sid"
  local ps="$BASEDIR/05_plasrisk_short/$sid"
  local pl="$BASEDIR/05_plasrisk_long/$sid"
  local circ="$BASEDIR/06_circularity/${sid}.circular.tsv"
  local paf="$BASEDIR/07_align/${sid}.paf"
  {
    echo "[truth] $sid start"
    # 1. long-read truth assembly
    if [[ ! -s "$long_fa" ]]; then
      download_url "$long_gz" "$url_long" "$BASEDIR/01_truth_longread" || exit 11
      gunzip -c "$long_gz" > "$long_fa" || exit 11
    fi
    # 2. matched Illumina reads
    mkdir -p "$rdir"
    download_url "$r1gz" "$r1" "$rdir" || exit 12
    download_url "$r2gz" "$r2" "$rdir" || exit 12
    # 3. short-read assembly (MEGAHIT)
    # MEGAHIT's --tmp-dir parent must already exist (tempfile.mkdtemp does
    # not create it). A failed attempt also leaves an output dir that
    # MEGAHIT refuses to reuse, so remove partial dirs before re-running.
    if [[ ! -s "$short_fa" ]]; then
      rm -rf "$mh" "$BASEDIR/03_megahit/tmp_$sid"
      mkdir -p "$BASEDIR/03_megahit/tmp_$sid"
      megahit -1 "$r1gz" -2 "$r2gz" -o "$mh" \
        --min-contig-len "$MIN_CONTIG" -t "$THREADS" \
        --tmp-dir "$BASEDIR/03_megahit/tmp_$sid" || exit 13
    fi
    # 4. geNomad + plasrisk meta on short-read assembly
    if [[ -z "$(find "$gs" -type f -name '*_plasmid_summary.tsv' 2>/dev/null | head -1)" ]]; then
      genomad end-to-end --cleanup --threads "$THREADS" "$short_fa" "$gs" "$GENOMAD_DB" || exit 14
    fi
    if [[ ! -s "$ps/meta_bin_scores.tsv" ]]; then
      plasrisk meta "$short_fa" --genomad-bin "$gs" $GE_ARG \
        --min-contig "$MIN_CONTIG" --complete-length 10000 \
        --threads "$THREADS" --strict-annotation -o "$ps" || exit 15
    fi
    # 5. geNomad + plasrisk meta on long-read truth assembly
    if [[ -z "$(find "$gl" -type f -name '*_plasmid_summary.tsv' 2>/dev/null | head -1)" ]]; then
      genomad end-to-end --cleanup --threads "$THREADS" "$long_fa" "$gl" "$GENOMAD_DB" || exit 16
    fi
    if [[ ! -s "$pl/meta_bin_scores.tsv" ]]; then
      plasrisk meta "$long_fa" --genomad-bin "$gl" $GE_ARG \
        --min-contig "$MIN_CONTIG" --complete-length 10000 \
        --threads "$THREADS" --strict-annotation -o "$pl" || exit 17
    fi
    # 6. circularity of long-read molecules (header flag + terminal overlap)
    if [[ ! -s "$circ" ]]; then
      python3 - "$long_fa" "$circ" "$CIRC_MINBP" "$AUTHOR_TABLES" <<'PYEOF'
import os, re, subprocess, sys
from collections import defaultdict
fasta, out, minbp, author_dir = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
def read_fasta(path):
    sid, buf = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if sid is not None:
                    yield sid, "".join(buf)
                sid = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if sid is not None:
            yield sid, "".join(buf)
hdr_re = re.compile(r"circular|_circ\b|_c$|circular=true", re.I)
# author Table C labels (optional): contig name -> circular
author = {}
for fn in os.listdir(author_dir):
    if fn.endswith(".done"):
        continue
    p = os.path.join(author_dir, fn)
    try:
        with open(p, encoding="utf-8", errors="replace") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ci = next((i for i, c in enumerate(head)
                       if re.search(r"circ|type|status", c, re.I)), None)
            ni = next((i for i, c in enumerate(head)
                       if re.search(r"contig|name|id", c, re.I)), 0)
            if ci is None:
                continue
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(ci, ni):
                    author[f[ni]] = f[ci].lower()
    except Exception:
        continue
# terminal-overlap: self-self minimap2 ava-ont, look for end-to-end overlaps
overlap = defaultdict(bool)
try:
    proc = subprocess.run(
        ["minimap2", "-x", "ava-ont", "-t", "4", fasta, fasta],
        capture_output=True, text=True, check=False)
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 10 or f[0] != f[5]:
            continue
        qs, qe, ts, te, L = int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[6])
        nm = int(f[9])
        idy = nm / max(1, (qe - qs) + nm)
        if idy < 0.90:
            continue
        q_wraps = (qs < 500 and te > L - 500) or (qe > L - 500 and ts < 500)
        if q_wraps and (qe - qs) >= minbp:
            overlap[f[0]] = True
except Exception as exc:
    print("[warn] minimap2 self-align failed: %s" % exc, file=sys.stderr)
with open(out, "w", encoding="utf-8") as oh:
    oh.write("seq_id\tlength_bp\theader_circular\toverlap_circular\tauthor_label\tcircular\n")
    for name, seq in read_fasta(fasta):
        h = bool(hdr_re.search(name))
        o = overlap.get(name, False)
        a = author.get(name, "")
        a_circ = "circ" in a
        circ = h or o or (a_circ if a else False)
        if len(seq) >= 1000:
            oh.write("%s\t%d\t%d\t%d\t%s\t%d\n"
                     % (name, len(seq), int(h), int(o), a or "NA", int(circ)))
PYEOF
    fi
    # 7. recovery alignment: short contigs (query) vs truth (target)
    if [[ ! -s "$paf" ]]; then
      minimap2 -x asm10 -c --secondary=no -t "$THREADS" \
        "$long_fa" "$short_fa" > "$paf" 2>/dev/null || exit 18
    fi
    echo "[truth] $sid DONE"
  } > "$lg" 2>&1
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[truth] $sid FAILED (exit $rc)"
    echo -e "${sid}\texit_${rc}" >> "${FAIL_LOG}.tmp"
    tail -n 20 "$lg"
  else
    echo "[truth] $sid DONE"
  fi
}

slot() { while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do sleep 3; done; }

# join the two manifests by id (id is the long-read alias)
mapfile -t ALL_ROWS < <(python3 - "$TRUTH_MANIFEST" "$SHORT_MANIFEST" <<'PYEOF'
import sys
longm, shortm = sys.argv[1], sys.argv[2]
long = {}
with open(longm, encoding="utf-8") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[0].startswith("GS_"):
            long[f[0]] = f[3]
short = {}
with open(shortm, encoding="utf-8") as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5:
            short[f[0]] = (f[3], f[4])
for sid, url_long in sorted(long.items()):
    if sid not in short:
        print("[warn] no matched Illumina for %s" % sid, file=sys.stderr)
        continue
    print("%s\t%s\t%s\t%s" % (sid, url_long, short[sid][0], short[sid][1]))
PYEOF
)

for row in "${ALL_ROWS[@]}"; do
  IFS=$'\t' read -r sid url_long r1 r2 <<< "$row"
  slot
  process_one "$sid" "$url_long" "$r1" "$r2" &
done
wait
mv "${FAIL_LOG}.tmp" "$FAIL_LOG" 2>/dev/null || true
echo "failed samples: $(wc -l < "$FAIL_LOG" 2>/dev/null || echo 0)"

# ---------------- funnel aggregation --------------------------------------------
python3 - "$BASEDIR" "$OVL_IDY" "$OVL_COV_STRICT" "$OVL_COV_RELAXED" "$MIN_CONTIG" <<'PYEOF'
import os, sys, glob
import pandas as pd

base, idy_thr, cov_strict, cov_relaxed, minc = sys.argv[1:6]
idy_thr, cov_strict, cov_relaxed, minc = float(idy_thr), float(cov_strict), float(cov_relaxed), int(minc)

def parse_paf(path):
    """Return dict truth_id -> covered_bases (identity-filtered)."""
    cov = {}
    best = {}
    if not os.path.exists(path):
        return cov, best
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, t, tlen = f[0], f[5], int(f[6])
            ts, te, nm = int(f[7]), int(f[8]), int(f[9])
            blen = int(f[10])
            idy = (blen - nm) / blen if blen else 0.0
            if idy < idy_thr:
                continue
            a = cov.setdefault(t, [0] * tlen)
            for i in range(ts, min(te, tlen)):
                a[i] = 1
            if t not in best or idy > best[t][0]:
                best[t] = (idy, q, (te - ts) / tlen)
    cov = {t: (sum(v), len(v)) for t, v in cov.items()}
    return cov, best

def called_contig_set(gdir):
    """seq names listed in geNomad plasmid_summary coordinates."""
    s = set()
    for p in glob.glob(os.path.join(gdir, "**", "*_plasmid_summary.tsv"),
                       recursive=True):
        try:
            d = pd.read_csv(p, sep="\t", dtype=str)
        except Exception:
            continue
        for _, r in d.iterrows():
            for col in ("seq_name", "coordinates"):
                if col in d.columns and isinstance(r.get(col), str):
                    for token in r[col].replace(",", " ").split():
                        name = token.split(":")[0].split("-")[0]
                        if name:
                            s.add(name)
            if "seq_name" in d.columns and isinstance(r.get("seq_name"), str):
                s.add(r["seq_name"])
    return s

def contig_grade(pdir):
    """contig_id -> grade (fallback to bins); only complete grades."""
    m = {}
    cp = os.path.join(pdir, "meta_contig_scores.tsv")
    if os.path.exists(cp):
        d = pd.read_csv(cp, sep="\t", dtype=str)
        gcol = "grade" if "grade" in d.columns else None
        scol = next((c for c in ("seq_id", "seq_name") if c in d.columns), None)
        if gcol and scol:
            for _, r in d.iterrows():
                if isinstance(r[gcol], str) and r[gcol]:
                    m[r[scol]] = r[gcol]
    return m

samples = sorted(os.path.basename(p) for p in
                 glob.glob(os.path.join(base, "01_truth_longread", "*.fna")))
funnel_rows, rec_rows, long_grade_rows = [], [], []

for sid in [s[:-4] for s in samples]:
    long_fa = os.path.join(base, "01_truth_longread", sid + ".fna")
    circ_p = os.path.join(base, "06_circularity", sid + ".circular.tsv")
    paf_p = os.path.join(base, "07_align", sid + ".paf")
    gs = os.path.join(base, "04_genomad_short", sid)
    ps = os.path.join(base, "05_plasrisk_short", sid)
    pl = os.path.join(base, "05_plasrisk_long", sid)
    if not (os.path.exists(circ_p) and os.path.exists(paf_p)):
        continue
    circ = pd.read_csv(circ_p, sep="\t")
    truth = circ[circ["circular"] == 1].copy()
    truth = truth[truth["length_bp"] >= minc]
    cov, best = parse_paf(paf_p)
    called = called_contig_set(gs)
    gmap = contig_grade(ps)
    # rebuild grade only for called contigs (bin-level rows handled below)
    for i, r in truth.iterrows():
        t = r["seq_id"]
        c, L = cov.get(t, (0, r["length_bp"]))
        frac = c / L if L else 0.0
        idy, qcontig, qfrac = (best.get(t, (None, None, None)) + (None,)*3)[:3]
        assembled = frac >= cov_strict
        relaxed = frac >= cov_relaxed
        is_called = bool(qcontig) and (qcontig in called)
        grade = gmap.get(qcontig, "") if qcontig else ""
        rec_rows.append({
            "sample": sid, "truth_id": t, "length_bp": L,
            "covered_frac": round(frac, 4), "best_idy": round(idy, 4) if idy else "",
            "best_query": qcontig or "", "best_query_frac": round(qfrac, 4) if qfrac else "",
            "recovered_strict": int(assembled), "recovered_relaxed": int(relaxed),
            "genomad_called": int(is_called), "plasrisk_grade": grade})
    n = len(truth)
    bins = [0, 1500, 2000, 3000, 5000, 10000, 10**9]
    labels = ["1-1.5k", "1.5-2k", "2-3k", "3-5k", "5-10k", ">=10k"]
    fr = {"sample": sid, "n_truth_circular": n}
    for lo, hi, lab in zip(bins[:-1], bins[1:], labels):
        sub = truth[(truth["length_bp"] >= lo) & (truth["length_bp"] < hi)]
        rec = pd.DataFrame(rec_rows)
        rec = rec[rec["sample"] == sid]
        m = rec.set_index("truth_id")
        ss = sub[sub["seq_id"].isin(m.index)]
        if len(ss):
            mm = m.loc[ss["seq_id"]]
            fr["n_%s" % lab] = len(ss)
            fr["rec_%s" % lab] = int(mm["recovered_strict"].sum())
    rdf = pd.DataFrame(rec_rows)
    if len(rdf):
        rdf = rdf[rdf["sample"] == sid]
        fr["n_recovered"] = int(rdf["recovered_strict"].sum())
        fr["n_recovered_relaxed"] = int(rdf["recovered_relaxed"].sum())
        fr["n_called"] = int(rdf["genomad_called"].sum())
        fr["n_graded"] = int((rdf["plasrisk_grade"] != "").sum())
        for g in ("A", "B", "C", "D", "E"):
            fr["n_grade_%s" % g] = int((rdf["plasrisk_grade"] == g).sum())
    funnel_rows.append(fr)
    # direct grades of long-read truth molecules
    bp = os.path.join(pl, "meta_bin_scores.tsv")
    if os.path.exists(bp):
        d = pd.read_csv(bp, sep="\t")
        if not d.empty:
            d["sample"] = sid
            long_grade_rows.append(d)

funnel = pd.DataFrame(funnel_rows)
funnel.to_csv(os.path.join(base, "truth_funnel.tsv"), sep="\t", index=False)
pd.DataFrame(rec_rows).to_csv(os.path.join(base, "truth_recovery_contigs.tsv"),
                              sep="\t", index=False)
if long_grade_rows:
    lg = pd.concat(long_grade_rows, ignore_index=True)
    keep = ["sample", "bin_id", "bin_origin", "total_length_bp", "confidence",
            "n_ARG", "n_VF", "n_BM", "high_risk_genes", "plasmid_score",
            "S_lite", "S_frag4", "S_cargo3", "grade", "grade_label"]
    lg[[c for c in keep if c in lg.columns]].to_csv(
        os.path.join(base, "truth_long_grades.tsv"), sep="\t", index=False)

# report
lines = []
if len(funnel):
    n_truth = int(funnel["n_truth_circular"].sum())
    n_rec = int(funnel["n_recovered"].sum())
    n_recx = int(funnel["n_recovered_relaxed"].sum())
    n_call = int(funnel["n_called"].sum())
    n_grad = int(funnel["n_graded"].sum())
    lines.append("samples: %d" % len(funnel))
    lines.append("truth circular molecules (>= %d bp): %d" % (minc, n_truth))
    lines.append("recovered by MEGAHIT (idy>=%.2f, cov>=%.0f%%): %d (%.1f%%)"
                 % (idy_thr, cov_strict*100, n_rec, 100*n_rec/max(1, n_truth)))
    lines.append("recovered (relaxed cov>=%.0f%%): %d (%.1f%%)"
                 % (cov_relaxed*100, n_recx, 100*n_recx/max(1, n_truth)))
    lines.append("recovered AND geNomad plasmid-called: %d (%.1f%% of truth)"
                 % (n_call, 100*n_call/max(1, n_truth)))
    lines.append("recovered AND PlasRisk graded: %d (%.1f%% of truth)"
                 % (n_grad, 100*n_grad/max(1, n_truth)))
    if long_grade_rows:
        lg = pd.concat(long_grade_rows)
        pred = lg[lg["bin_origin"] == "predicted"] if "bin_origin" in lg.columns else lg
        g = pred["grade"].value_counts().to_dict() if "grade" in pred.columns else {}
        lines.append("direct PlasRisk grades on long-read predicted bins: "
                     + ", ".join("%s=%d" % (k, g.get(k, 0)) for k in ("A","B","C","D","E")))
        if "high_risk_genes" in pred.columns:
            lines.append("long-read predicted bins with high-risk ARG: %d"
                         % int((pred["high_risk_genes"].fillna("").astype(str) != "").sum()))
txt = "\n".join(lines)
print(txt)
with open(os.path.join(base, "truth_summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(txt + "\n")
PYEOF

echo "[truth] done -> $BASEDIR/truth_summary.txt"
