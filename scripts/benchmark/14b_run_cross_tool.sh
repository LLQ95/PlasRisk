#!/usr/bin/env bash
# 14b_run_cross_tool.sh
#
# Run-only part of the cross-tool comparison. Assumes the two conda
# environments and the databases are already installed via 14a_setup_envs.sh;
# no environment is created or modified here.
#
# P0-3 head-to-head on the 189 traveller diarrheagenic E. coli isolates
# (PRJNA698223): RFPlasmid 0.0.18 (Enterobacteriaceae model) and Plasmer
# are run on the same assemblies used by 06, and per-contig calls are joined
# with geNomad calls and PlasRisk grades.
#
# Usage:
#   bash 14b_run_cross_tool.sh results/traveler_isolates [threads] [rfp_jobs] [plm_jobs]
#
#   threads   threads handed to the tools (default 16; note RFPlasmid 0.0.18
#             hardcodes -t 16 for CheckM/jellyfish and lets diamond use all
#             cores, so keep rfp_jobs modest)
#   rfp_jobs  parallel RFPlasmid isolates (default 4)
#   plm_jobs  parallel Plasmer isolates (default 1; each kmer-db job needs
#             ~32 GB RAM, do not raise unless the machine has the memory)
#
# The join step runs in the CALLING environment and needs pandas, so run
# this script from the genomad conda env:
#   conda activate genomad
#   bash 14b_run_cross_tool.sh results/traveler_isolates 16 4 1
#
# Outputs (under <iso>/cross_tool/):
#   rfp_inputs/<acc>/<acc>.fasta      (one-symlink input dirs for RFPlasmid)
#   rfplasmid/<acc>/prediction.csv
#   plasmer/<acc>/results/<acc>.plasmer.predClass.tsv (+predProb.tsv)
#   logs/<acc>.rfp.log, logs/<acc>.plm.log
#   tool_run_status.tsv
#   tool_comparison_contigs.tsv       (per-contig four-tool join, all contigs)
#   tool_comparison_summary.tsv       (call counts per tool and grade set)
#   tool_comparison_agreement.tsv     (two-way agreement vs geNomad)
#   tool_comparison_critical.tsv      (the eight critical-ARG contigs)
#
# RFPlasmid quirks handled here (verified against rfplasmid 0.0.18 source):
#   * --input must be a directory; only *.fasta files are read (*.fna is
#     ignored), so each isolate gets its own directory with a .fasta symlink.
#   * --out must NOT pre-exist: if it exists RFPlasmid appends a
#     _YYYYMMDD_HHMMSS suffix and writes there instead.
#   * prediction.csv columns: prediction (c/p), votes chromosomal,
#     votes plasmid, contigID (the original full FASTA header).
# Plasmer quirks handled here (verified against nekokoe/Plasmer main):
#   * launcher is `Plasmer -g fa -p prefix -d db -t t -m min -l length -o out`
#   * -l 0 is required so every contig is predicted (default 500000 auto-calls
#     long contigs chromosome without modelling); -m 500 leaves short contigs
#     labelled shorter_than_500, counted as unclassified.
#   * success is judged by the presence of predClass.tsv, not the exit code
#     (the launcher exits 1 when the optional kraken DB is absent).

set -uo pipefail

# English locale for all log output (timestamps, tool messages).
export LC_ALL=C

ISO_DIR="${1:-results/traveler_isolates}"
THREADS="${2:-16}"
RFP_JOBS="${3:-4}"
PLM_JOBS="${4:-1}"
ASM_DIR="$ISO_DIR/assemblies"
OUTDIR="$ISO_DIR/cross_tool"
ENV_RFP="rfplasmid_cmp"
ENV_PLM="plasmer_cmp"
PLASMER_DB="${PLASMER_DB:-$OUTDIR/db/plasmer_db}"

mkdir -p "$OUTDIR/rfp_inputs" "$OUTDIR/logs"

# ------------------------------------------------------------- prechecks
[[ -d "$ASM_DIR" ]] || { echo "[error] $ASM_DIR not found (run 06 first)"; exit 1; }
command -v conda >/dev/null 2>&1 || { echo "[error] conda not in PATH"; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env list | awk '{print $1}' | grep -qx "$ENV_RFP" || {
  echo "[error] env $ENV_RFP missing; run 14a_setup_envs.sh first"; exit 2; }
conda env list | awk '{print $1}' | grep -qx "$ENV_PLM" || {
  echo "[error] env $ENV_PLM missing; run 14a_setup_envs.sh first"; exit 2; }
[[ -s "$PLASMER_DB/plasmerMainDB.list" ]] || {
  echo "[error] $PLASMER_DB/plasmerMainDB.list missing; run 14a_setup_envs.sh"
  echo "        (or export PLASMER_DB=/path/to/extracted/db)"; exit 3; }
[[ -s "$PLASMER_DB/customizedKraken2DB.list" ]] || {
  echo "[error] $PLASMER_DB/customizedKraken2DB.list missing; run 14a_setup_envs.sh"; exit 3; }
python3 -c "import pandas" 2>/dev/null || {
  echo "[error] the calling env needs pandas (run from the genomad env)"; exit 4; }

# RFPlasmid model DB readiness (both .dmnd files inside the package)
RFP_READY=0
conda activate "$ENV_RFP"
RFP_DIR="$(python3 -c 'import RFPlasmid, os; print(os.path.dirname(RFPlasmid.__file__))' 2>/dev/null)"
if [[ -n "$RFP_DIR" && -s "$RFP_DIR/plasmiddb_cge.dmnd" \
      && -s "$RFP_DIR/plasmiddb_total.dmnd" ]]; then
  RFP_READY=1
  echo "[cmp] RFPlasmid DB ready: $RFP_DIR"
else
  echo "[error] RFPlasmid model DB not initialised; re-run 14a_setup_envs.sh"
fi
conda deactivate
# Plasmer launcher must be on PATH inside the plasmer env
conda activate "$ENV_PLM"
command -v Plasmer >/dev/null 2>&1 || {
  echo "[error] Plasmer launcher not on PATH inside $ENV_PLM";
  conda deactivate; exit 5; }
echo "[cmp] Plasmer launcher ready: $(command -v Plasmer)"
conda deactivate
[[ "$RFP_READY" == 1 ]] || exit 5

n_assemblies=$(find "$ASM_DIR" -maxdepth 1 -name '*.fna' | wc -l)
echo "[cmp] $n_assemblies assemblies, threads=$THREADS rfp_jobs=$RFP_JOBS plm_jobs=$PLM_JOBS"

# ------------------------------------------------------- RFPlasmid phase
run_rfp() {
  local acc="$1"
  local fa="$ASM_DIR/${acc}.fna"
  local idir="$OUTDIR/rfp_inputs/$acc"
  local rdir="$OUTDIR/rfplasmid/$acc"
  local lg="$OUTDIR/logs/${acc}.rfp.log"
  # already done in the canonical output dir?
  if [[ -s "$rdir/prediction.csv" ]]; then
    echo "[rfp] $acc already done"; return 0
  fi
  # already done in a timestamped dir (created by a previous run)?
  if find "$OUTDIR/rfplasmid" -maxdepth 2 -type f -name prediction.csv \
       \( -path "*/${acc}/prediction.csv" \
          -o -path "*/${acc}_[0-9]*_[0-9]*/prediction.csv" \) 2>/dev/null | grep -q .; then
    echo "[rfp] $acc already done (timestamped dir)"; return 0
  fi
  {
    echo "[rfp] $acc start $(date)"
    # per-isolate input directory containing exactly one *.fasta symlink
    mkdir -p "$idir"
    ln -sfn "$(readlink -f "$fa")" "$idir/${acc}.fasta"
    # RFPlasmid refuses to reuse an existing --out and appends a suffix;
    # remove only an empty/stale dir that has no prediction
    if [[ -d "$rdir" && -z "$(find "$rdir" -type f -name 'prediction.csv' 2>/dev/null | head -1)" ]]; then
      rm -rf "$rdir"
    fi
    conda activate "$ENV_RFP"
    rfplasmid --species Enterobacteriaceae --input "$idir" \
      --out "$rdir" --threads "$THREADS"
    rc=$?
    conda deactivate
    echo "[rfp] $acc exit=$rc $(date)"
  } > "$lg" 2>&1
  if [[ -s "$rdir/prediction.csv" ]] || \
     find "$OUTDIR/rfplasmid" -maxdepth 2 -type f -name prediction.csv \
       \( -path "*/${acc}/prediction.csv" \
          -o -path "*/${acc}_[0-9]*_[0-9]*/prediction.csv" \) 2>/dev/null | grep -q .; then
    echo "[rfp] $acc DONE"
  else
    echo "[rfp] $acc FAILED (see $lg)"
  fi
}

# --------------------------------------------------------- Plasmer phase
run_plm() {
  local acc="$1"
  local fa="$ASM_DIR/${acc}.fna"
  local pdir="$OUTDIR/plasmer/$acc"
  local lg="$OUTDIR/logs/${acc}.plm.log"
  if [[ -s "$pdir/results/${acc}.plasmer.predClass.tsv" ]]; then
    echo "[plm] $acc already done"; return 0
  fi
  {
    echo "[plm] $acc start $(date)"
    conda activate "$ENV_PLM"
    # -l 0: predict every contig; -m 500: below 500 bp -> shorter_than_500
    Plasmer -g "$fa" -p "$acc" -d "$PLASMER_DB" \
      -t "$THREADS" -m 500 -l 0 -o "$pdir"
    rc=$?
    conda deactivate
    # launcher may exit 1 after finishing the taxonomy step; trust outputs
    echo "[plm] $acc exit=$rc $(date)"
  } > "$lg" 2>&1
  if [[ -s "$pdir/results/${acc}.plasmer.predClass.tsv" ]]; then
    echo "[plm] $acc DONE"
  else
    echo "[plm] $acc FAILED (see $lg)"
  fi
}

slot() {
  local max="$1"
  while [[ "$(jobs -rp | wc -l)" -ge "$max" ]]; do sleep 3; done
}

n=0
for fa in "$ASM_DIR"/*.fna; do
  acc="$(basename "$fa" .fna)"
  slot "$RFP_JOBS"
  run_rfp "$acc" &
  n=$((n+1))
done
wait
echo "[cmp] RFPlasmid phase processed $n isolates"

n=0
for fa in "$ASM_DIR"/*.fna; do
  acc="$(basename "$fa" .fna)"
  slot "$PLM_JOBS"
  run_plm "$acc" &
  n=$((n+1))
done
wait
echo "[cmp] Plasmer phase processed $n isolates"

# ---------------- join everything ----------------------------------------------
python3 - "$ISO_DIR" "$OUTDIR" <<'PYEOF'
import glob, os, re, sys
import pandas as pd

iso, out = sys.argv[1], sys.argv[2]
asm = os.path.join(iso, "assemblies")


def norm_id(x):
    if not isinstance(x, str):
        return x
    return x.strip().lstrip(">").split()[0]


# ---- universe: every contig in every assembly
univ_rows = []
for fp in sorted(glob.glob(os.path.join(asm, "*.fna"))):
    acc = os.path.basename(fp)[:-4]
    with open(fp, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                univ_rows.append((acc, norm_id(line[1:])))
univ = pd.DataFrame(univ_rows, columns=["accession", "seq_id"])
print("[join] universe contigs:", len(univ))

# ---- PlasRisk grades (geNomad-called plasmid contigs scored by plasrisk)
sc_p = os.path.join(iso, "benchmark_isolate_scores.tsv")
sc = pd.read_csv(sc_p, sep="\t", dtype={"seq_id": str})
sc["seq_id"] = sc["seq_id"].map(norm_id)
keep = [c for c in ("accession", "seq_id", "S_norm", "grade",
                    "high_risk_genes", "n_ARG", "length_bp") if c in sc.columns]
sc = sc[keep].copy()
sc["plasrisk_scored"] = "yes"

# ---- geNomad per-contig summaries (only plasmid-called sequences are in them)
# Accession is taken from the summary file name (<acc>_plasmid_summary.tsv),
# because contig ids such as NODE_1 repeat across assemblies.
ge_rows = []
for p in glob.glob(os.path.join(iso, "genomad", "**", "*_plasmid_summary.tsv"),
                  recursive=True):
    d = pd.read_csv(p, sep="\t", dtype=str)
    if "seq_name" not in d.columns:
        continue
    acc = os.path.basename(p).replace("_plasmid_summary.tsv", "")
    d["accession"] = acc
    d["seq_id"] = d["seq_name"].map(norm_id)
    cols = ["accession", "seq_id"] + [
        c for c in ("plasmid_score", "length", "topology") if c in d.columns]
    ge_rows.append(d[cols])
ge = (pd.concat(ge_rows, ignore_index=True).drop_duplicates(
    ["accession", "seq_id"]) if ge_rows else
    pd.DataFrame(columns=["accession", "seq_id"]))
ge = ge.rename(columns={"plasmid_score": "genomad_plasmid_score"})
ge["genomad_call"] = "plasmid"

# ---- RFPlasmid: prediction.csv has columns
#       index(<genome>_<n>), prediction(c/p), votes chromosomal,
#       votes plasmid, contigID (original full header)
rfp_rows = []
for fp in glob.glob(os.path.join(out, "rfplasmid", "**", "prediction.csv"),
                   recursive=True):
    acc = os.path.basename(os.path.dirname(fp))
    acc = re.sub(r"_\d{8}_\d{6}$", "", acc)  # strip _YYYYMMDD_HHMMSS suffix
    d = pd.read_csv(fp, index_col=0)
    d = d.rename(columns={"prediction": "rfplasmid_call",
                          "votes plasmid": "rfplasmid_plasmid_prob"})
    d["seq_id"] = d["contigID"].map(norm_id)
    d["accession"] = acc
    rfp_rows.append(d[["accession", "seq_id", "rfplasmid_call",
                       "rfplasmid_plasmid_prob"]])
rfp = (pd.concat(rfp_rows, ignore_index=True).drop_duplicates(
    ["accession", "seq_id"]) if rfp_rows else
    pd.DataFrame(columns=["accession", "seq_id", "rfplasmid_call",
                          "rfplasmid_plasmid_prob"]))
rfp["rfplasmid_call"] = rfp["rfplasmid_call"].map(
    {"p": "plasmid", "c": "chromosome"}).fillna(rfp["rfplasmid_call"])

# ---- Plasmer: predClass.tsv is headerless (id, class); predProb.tsv has
#      Contig/chromosome/plasmid
plm_rows = []
for d_acc in sorted(glob.glob(os.path.join(out, "plasmer", "*"))):
    acc = os.path.basename(d_acc)
    fcls = os.path.join(d_acc, "results", acc + ".plasmer.predClass.tsv")
    if not os.path.exists(fcls):
        continue
    c = pd.read_csv(fcls, sep="\t", header=None,
                    names=["seq_id", "plasmer_call"], dtype=str)
    c["seq_id"] = c["seq_id"].map(norm_id)
    c["accession"] = acc
    fprob = os.path.join(d_acc, "results", acc + ".plasmer.predProb.tsv")
    if os.path.exists(fprob):
        pr = pd.read_csv(fprob, sep="\t", dtype=str)
        # Three real header layouts exist across Plasmer releases:
        #   "Contig chromosome plasmid"       documented; 3 names
        #   "<tab>chromosome plasmid"          leading tab; ids -> Unnamed: 0
        #   "chromosome plasmid"               2 names for 3 fields; ids -> index
        if "Contig" in pr.columns:
            pr = pr.rename(columns={"Contig": "seq_id"})
        elif "Unnamed: 0" in pr.columns:
            pr = pr.rename(columns={"Unnamed: 0": "seq_id"})
        else:
            pr = pr.reset_index(names="seq_id")
        pr = pr.rename(columns={"plasmid": "plasmer_plasmid_prob"})
        pr["seq_id"] = pr["seq_id"].map(norm_id)
        pr = pr.drop_duplicates("seq_id")
        c = c.merge(pr[["seq_id", "plasmer_plasmid_prob"]], on="seq_id",
                    how="left")
    plm_rows.append(c)
plm = (pd.concat(plm_rows, ignore_index=True).drop_duplicates(
    ["accession", "seq_id"]) if plm_rows else
    pd.DataFrame(columns=["accession", "seq_id", "plasmer_call",
                          "plasmer_plasmid_prob"]))


def norm_class(x):
    if not isinstance(x, str):
        return ""
    x = x.lower()
    if x.startswith("shorter"):
        return "unclassified"
    if "chrom" in x:
        return "chromosome"
    if "plas" in x or x == "p":
        return "plasmid"
    return x


plm["plasmer_call"] = plm["plasmer_call"].map(norm_class)

# ---- merge
m = univ.merge(sc, on=["accession", "seq_id"], how="left")
m = m.merge(ge, on=["accession", "seq_id"], how="left")
m = m.merge(rfp, on=["accession", "seq_id"], how="left")
m = m.merge(plm, on=["accession", "seq_id"], how="left")
m["genomad_call"] = m["genomad_call"].fillna("chromosome")
m.to_csv(os.path.join(out, "tool_comparison_contigs.tsv"),
         sep="\t", index=False)

# ---- run status
status = m.groupby("accession").agg(
    n_contigs=("seq_id", "size"),
    rfp_classified=("rfplasmid_call", lambda s: int(s.isin(
        ["plasmid", "chromosome"]).sum())),
    plm_classified=("plasmer_call", lambda s: int(s.isin(
        ["plasmid", "chromosome"]).sum())),
    plm_unclassified=("plasmer_call", lambda s: int((s == "unclassified").sum())),
).reset_index()
status["rfp_ok"] = status["rfp_classified"] > 0
status["plm_ok"] = status["plm_classified"] > 0
status.to_csv(os.path.join(out, "tool_run_status.tsv"), sep="\t", index=False)
print("[join] RFPlasmid ok isolates:", int(status["rfp_ok"].sum()),
      "/", len(status))
print("[join] Plasmer ok isolates:", int(status["plm_ok"].sum()),
      "/", len(status))

# ---- summaries
summary = []
for tool, col in [("RFPlasmid", "rfplasmid_call"),
                  ("Plasmer", "plasmer_call")]:
    def block(name, df):
        v = df[col]
        summary.append({
            "tool": tool, "set": name, "n": len(df),
            "called_plasmid": int((v == "plasmid").sum()),
            "called_chromosome": int((v == "chromosome").sum()),
            "unclassified": int((v == "unclassified").sum()),
            "missing": int(v.isna().sum())})
    block("all_assembly_contigs", m)
    block("genomad_plasmid", m[m["genomad_call"] == "plasmid"])
    for g in ["A", "B", "C", "D", "E"]:
        block("plasrisk_grade_%s" % g, m[m["grade"] == g])
    block("plasrisk_grade_A_or_B", m[m["grade"].isin(["A", "B"])])
    block("plasrisk_grade_C_or_above", m[m["grade"].isin(["A", "B", "C"])])
pd.DataFrame(summary).to_csv(
    os.path.join(out, "tool_comparison_summary.tsv"), sep="\t", index=False)

# ---- two-way agreement against geNomad on the full contig universe
agree = []
for tool, col in [("RFPlasmid", "rfplasmid_call"),
                  ("Plasmer", "plasmer_call")]:
    v = m[col]
    g = m["genomad_call"]
    tp = int(((v == "plasmid") & (g == "plasmid")).sum())
    fp = int(((v == "plasmid") & (g == "chromosome")).sum())
    fn = int(((v == "chromosome") & (g == "plasmid")).sum())
    tn = int(((v == "chromosome") & (g == "chromosome")).sum())
    agree.append({"tool": tool, "reference": "geNomad",
                  "both_plasmid": tp, "both_chromosome": tn,
                  "tool_only_plasmid": fp, "geNomad_only_plasmid": fn,
                  "total_classified_by_tool": int(v.isin(
                      ["plasmid", "chromosome"]).sum())})
pd.DataFrame(agree).to_csv(
    os.path.join(out, "tool_comparison_agreement.tsv"), sep="\t", index=False)

# ---- the eight critical-ARG contigs
crit_p = os.path.join(iso, "critical_gene_hits.tsv")
if os.path.exists(crit_p):
    crit = pd.read_csv(crit_p, sep="\t", dtype=str)
    ccol = next((c for c in ("seq_id", "contig", "contig_id")
                 if c in crit.columns), None)
    if ccol:
        crit["seq_id"] = crit[ccol].map(norm_id)
        idx = m.set_index(["accession", "seq_id"])
        cols = ["genomad_call", "genomad_plasmid_score", "rfplasmid_call",
                "rfplasmid_plasmid_prob", "plasmer_call",
                "plasmer_plasmid_prob", "grade", "S_norm"]
        for c in cols:
            if c in m.columns:
                crit[c] = [
                    (idx.loc[(a, s), c] if (a, s) in idx.index else "ABSENT")
                    for a, s in zip(crit["accession"], crit["seq_id"])]
        crit.to_csv(os.path.join(out, "tool_comparison_critical.tsv"),
                    sep="\t", index=False)

print(pd.DataFrame(summary).to_string(index=False))
print("written:", out)
PYEOF

echo "[cmp] done; see $OUTDIR/tool_comparison_*.tsv"
