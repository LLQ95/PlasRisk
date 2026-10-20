#!/usr/bin/env bash
# 14_rfplasmid_plasmer.sh
#
# P0-3 head-to-head on the 189 traveller diarrheagenic E. coli
# isolates (PRJNA698223): RFPlasmid and Plasmer are run on the same
# downloaded assemblies used in 06_traveler_isolates_benchmark.sh, and
# their per-contig calls are joined with geNomad calls and PlasRisk
# grades. RFPlasmid is constrained to the Enterobacteriaceae model,
# matching the target organisms.
#
# Inputs (already produced by 06):
#   results/traveler_isolates/assemblies/<accession>.fna
#   results/traveler_isolates/genomad/<acc>/.../*_plasmid_summary.tsv
#   results/traveler_isolates/benchmark_isolate_scores.tsv
#   results/traveler_isolates/critical_gene_hits.tsv
#
# Usage:
#   bash 14_rfplasmid_plasmer.sh results/traveler_isolates 8 8
#
# Note: the script installs two isolated conda environments on first
# run. Plasmer's database URL is taken from its GitHub README; if the
# hardcoded URL no longer resolves, the script prints the exact manual
# step and stops (it does not guess silently).
#
# Outputs:
#   outdir/rfplasmid/<acc>/, plasmer/<acc>/
#   tool_comparison_contigs.tsv   (per-contig four-tool join)
#   tool_comparison_summary.tsv   (agreement statistics)
#   tool_comparison_critical.tsv  (the eight critical-ARG contigs)

set -uo pipefail

ISO_DIR="${1:-results/traveler_isolates}"
THREADS="${2:-8}"
JOBS="${3:-8}"
ASM_DIR="$ISO_DIR/assemblies"
OUTDIR="$ISO_DIR/cross_tool"
mkdir -p "$OUTDIR"/{rfplasmid,plasmer,logs}

[[ -d "$ASM_DIR" ]] || { echo "[error] $ASM_DIR not found (run 06 first)"; exit 1; }
command -v conda >/dev/null 2>&1 || { echo "[error] conda not in PATH"; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"

CONDA_BASE="$(conda info --base)"
ENV_RFP="rfplasmid_cmp"
ENV_PLM="plasmer_cmp"
PLASMER_REPO="$OUTDIR/Plasmer"
PLASMER_DB_URL="https://zenodo.org/records/10201491/files/PlasmerDB.tar.gz"

# ---------------- env setup (idempotent) ---------------------------------------
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_RFP"; then
  echo "[setup] creating $ENV_RFP"
  conda create -y -n "$ENV_RFP" -c conda-forge -c bioconda \
    rfplasmid python=3.10 || exit 20
fi
conda activate "$ENV_RFP"
if ! rfplasmid --help 2>/dev/null | grep -qi "predictions\|prediction"; then
  # model database initialisation (one-time)
  rfplasmid --prepareDB || rfplasmid --initialize || {
    echo "[error] RFPlasmid database initialisation failed; run 'rfplasmid --prepareDB' manually"; exit 21; }
fi
conda deactivate

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_PLM"; then
  echo "[setup] creating $ENV_PLM"
  conda create -y -n "$ENV_PLM" -c conda-forge -c bioconda \
    python=3.10 "diamond=2.0.15" jellyfish prodigal scikit-learn pandas numpy requests || exit 22
  git clone --depth 1 https://github.com/tiagolbiotech/Plasmer.git "$PLASMER_REPO" || exit 23
fi
if [[ ! -d "$PLASMER_REPO/PlasmerDB" ]]; then
  echo "[setup] downloading PlasmerDB"
  if curl -L --retry 3 -o "$OUTDIR/PlasmerDB.tar.gz" "$PLASMER_DB_URL"; then
    tar -xzf "$OUTDIR/PlasmerDB.tar.gz" -C "$PLASMER_REPO" || exit 24
  else
    echo "[error] PlasmerDB URL not reachable. Open $PLASMER_REPO/README.md"
    echo "       and place the database as $PLASMER_REPO/PlasmerDB, then rerun."
    exit 24
  fi
fi

run_one() {
  local acc="$1"
  local fa="$ASM_DIR/${acc}.fna"
  local lg="$OUTDIR/logs/${acc}.log"
  {
    echo "[cmp] $acc start"
    conda activate "$ENV_RFP"
    local rdir="$OUTDIR/rfplasmid/$acc"
    if [[ -z "$(find "$rdir" -type f -name '*prediction*' 2>/dev/null | head -1)" ]]; then
      mkdir -p "$rdir"
      rfplasmid --species Enterobacteriaceae --input "$fa" \
        --output "$rdir" --threads "$THREADS" 2>&1 || \
      rfplasmid -s Enterobacteriaceae -i "$fa" -o "$rdir" -j "$THREADS" 2>&1
    fi
    conda deactivate
    conda activate "$ENV_PLM"
    local pdir="$OUTDIR/plasmer/$acc"
    if [[ -z "$(find "$pdir" -type f -name '*.csv' 2>/dev/null | head -1)" ]]; then
      mkdir -p "$pdir"
      python "$PLASMER_REPO/Plasmer.py" -q "$fa" -o "$pdir" \
        -d "$PLASMER_REPO/PlasmerDB" -t "$THREADS" 2>&1
    fi
    conda deactivate
    echo "[cmp] $acc DONE"
  } > "$lg" 2>&1
  if [[ $? -ne 0 ]]; then
    echo "[cmp] $acc FAILED (see $lg)"
  else
    echo "[cmp] $acc DONE"
  fi
}

slot() { while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do sleep 3; done; }
n=0
for fa in "$ASM_DIR"/*.fna; do
  acc="$(basename "$fa" .fna)"
  slot
  run_one "$acc" &
  n=$((n+1))
done
wait
echo "processed $n isolates"

# ---------------- join everything ----------------------------------------------
python3 - "$ISO_DIR" "$OUTDIR" <<'PYEOF'
import glob, os, re, sys
import pandas as pd

iso, out = sys.argv[1], sys.argv[2]

# PlasRisk grades (geNomad-called plasmid contigs)
sc = pd.read_csv(os.path.join(iso, "benchmark_isolate_scores.tsv"), sep="\t")
seq_col = next(c for c in ("seq_id", "seq_name", "contig") if c in sc.columns)
sc = sc.rename(columns={seq_col: "seq_id"})

# geNomad per-contig score (called plasmids + every contig in summaries)
rows = []
for p in glob.glob(os.path.join(iso, "genomad", "**", "*_plasmid_summary.tsv"),
                   recursive=True):
    d = pd.read_csv(p, sep="\t")
    if "seq_name" in d.columns:
        d = d.rename(columns={"seq_name": "seq_id"})
        rows.append(d[["seq_id"] + [c for c in ("plasmid_score", "length", "topology")
                                   if c in d.columns]])
ge = pd.concat(rows, ignore_index=True).drop_duplicates("seq_id") if rows else \
     pd.DataFrame(columns=["seq_id"])
ge["genomad_call"] = "plasmid"

def find_col(cols, *cands):
    for cand in cands:
        for c in cols:
            if cand in c.lower():
                return c
    return None

def read_tool_dir(tool, pattern):
    rows = []
    for acc_dir in sorted(glob.glob(os.path.join(out, tool, "*"))):
        files = [f for f in glob.glob(os.path.join(acc_dir, "**", "*"),
                                      recursive=True)
                 if re.search(pattern, os.path.basename(f))]
        for f in files:
            try:
                d = pd.read_csv(f, sep=None, engine="python", dtype=str)
            except Exception:
                continue
            if d.empty:
                continue
            idc = find_col(d.columns, "seq_id", "contig", "name", "sequence")
            if idc is None:
                continue
            d = d.rename(columns={idc: "seq_id"})
            d["accession"] = os.path.basename(acc_dir)
            keep = ["accession", "seq_id"]
            pc = find_col(d.columns, "plasmid", "prediction", "label", "class")
            if pc:
                keep.append(pc)
                d = d.rename(columns={pc: tool + "_call"})
            sc_ = find_col(d.columns, "plasmid_prob", "plasmid probability",
                           "probability", "score")
            if sc_:
                keep.append(sc_)
                d = d.rename(columns={sc_: tool + "_plasmid_prob"})
            rows.append(d[[c for c in keep if c in d.columns]])
    if not rows:
        return pd.DataFrame(columns=["accession", "seq_id"])
    return pd.concat(rows, ignore_index=True).drop_duplicates(
        ["accession", "seq_id"])

rf = read_tool_dir("rfplasmid", r"prediction")
pl = read_tool_dir("plasmer", r"\.csv$")

# attach accession to PlasRisk table
base = sc[["accession", "seq_id", "S_norm", "grade", "high_risk_genes",
           "n_ARG", "length_bp"]].copy()
merged = base.merge(ge, on="seq_id", how="left")
merged = merged.merge(rf, on=["accession", "seq_id"], how="left")
merged = merged.merge(pl, on=["accession", "seq_id"], how="left")
merged.to_csv(os.path.join(out, "tool_comparison_contigs.tsv"),
              sep="\t", index=False)

def norm(x):
    if not isinstance(x, str):
        return ""
    x = x.lower()
    if "chrom" in x:
        return "chromosome"
    if "plas" in x or x == "p":
        return "plasmid"
    return x

summary = []
for tool, col in [("RFPlasmid", "rfplasmid_call"),
                  ("Plasmer", "plasmer_call")]:
    if col not in merged.columns:
        continue
    call = merged[col].map(norm)
    plas = call == "plasmid"
    chrom = call == "chromosome"
    n = int((plas | chrom).sum())
    n_p = int(plas.sum())
    agree = int((plas & merged["genomad_call"].eq("plasmid")).sum())
    for grade in ["A", "B", "C", "D", "E"]:
        sub = merged[merged["grade"] == grade]
        c = sub[col].map(norm)
        summary.append({"tool": tool, "set": "grade_%s" % grade,
                        "n": len(sub),
                        "called_plasmid": int((c == "plasmid").sum()),
                        "called_chromosome": int((c == "chromosome").sum())})
    summary.append({"tool": tool, "set": "all_geNomad_plasmid_contigs",
                    "n": n, "called_plasmid": n_p,
                    "called_chromosome": int(chrom.sum())})
    summary.append({"tool": tool, "set": "grade_A_or_B",
                    "n": int(merged["grade"].isin(["A", "B"]).sum()),
                    "called_plasmid": int((plas & merged["grade"].isin(["A", "B"])).sum()),
                    "called_chromosome": int((chrom & merged["grade"].isin(["A", "B"])).sum())})
pd.DataFrame(summary).to_csv(os.path.join(out, "tool_comparison_summary.tsv"),
                             sep="\t", index=False)

# the eight critical-ARG contigs
crit_p = os.path.join(iso, "critical_gene_hits.tsv")
if os.path.exists(crit_p):
    crit = pd.read_csv(crit_p, sep="\t", dtype=str)
    ccol = next((c for c in ("seq_id", "contig", "contig_id") if c in crit.columns), None)
    if ccol:
        idx = merged.set_index(["accession", "seq_id"])
        for tool in ("genomad_call", "rfplasmid_call", "plasmer_call", "grade"):
            if tool in merged.columns:
                crit[tool] = [
                    (idx.loc[(a, c), tool] if (a, c) in idx.index else "ABSENT")
                    for a, c in zip(crit["accession"], crit[ccol])]
        crit.to_csv(os.path.join(out, "tool_comparison_critical.tsv"),
                    sep="\t", index=False)

print(pd.DataFrame(summary).to_string(index=False))
print("written:", out)
PYEOF
