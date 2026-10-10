#!/usr/bin/env bash
# =============================================================================
# 01_pilot_fragment_grid.sh
# Realistic-layer fragmentation grid for PlasRisk (Tier B)
#   closed plasmids -> synthetic communities -> simulated Illumina reads
#   -> assembly -> geNomad plasmid bins -> plasrisk meta -> matched truth
#
# Run on the server under any conda env that has: art/MEGAHIT/metaSPAdes,
# genomad, bowtie2, jgi_summarize_bam_contig_depths, plasrisk (>=1.5.0).
# Edit the CONFIG block before running.
# =============================================================================
set -euo pipefail

# ---------------------------- CONFIG (edit me) -------------------------------
PLASMID_FASTA="/db/student/metagenome/pipdb/contigs/truth/closed_plasmids.fasta"
TRUTH_TABLE="/db/student/metagenome/pipdb/contigs/truth/closed_plasmids_truth.tsv"
# truth table columns: plasmid_id,length_bp,S_whole,grade,high_risk_genes
WORKDIR="/db/student/metagenome/pipdb/contigs/fragment_grid"
THREADS=16
READ_LENGTHS=(150)            # e.g. 100 150 250
DEPTHS=(10 30)                # per-community total depth split across plasmids
REPLICATES=3
PLASMIDS_PER_COMMUNITY=30
ASSEMBLER="megahit"           # megahit | metaspades
# -----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"; cd "${WORKDIR}"

# Split closed plasmids into one-file-per-plasmid for community mixing
mkdir -p refs
if [ ! -f refs/.split_done ]; then
  awk '/^>/{n=$1; sub(/^>/,"",n); f="refs/"n".fa"} {print > f}' "${PLASMID_FASTA}"
  touch refs/.split_done
fi
mapfile -t ALL_PLASMIDS < <(ls refs/*.fa)
TOTAL_N=${#ALL_PLASMIDS[@]}
echo "[grid] ${TOTAL_N} closed reference plasmids available"

GRID_LOG="${WORKDIR}/grid_manifest.tsv"
echo -e "rep\tread_len\tdepth\tcommunity\tdir" > "${GRID_LOG}"

for rep in $(seq 1 "${REPLICATES}"); do
  for rl in "${READ_LENGTHS[@]}"; do
    for dep in "${DEPTHS}"; do
      tag="rep${rep}_r${rl}_d${dep}"
      rdir="${WORKDIR}/runs/${tag}"
      mkdir -p "${rdir}"; cd "${rdir}"

      # Deterministic random community of PLASMIDS_PER_COMMUNITY plasmids with
      # log-uniform relative abundances
      comm="community_${tag}.fa"
      python3 - "${WORKDIR}/refs" "${comm}" "${PLASMIDS_PER_COMMUNITY}" \
              "$((1000+rep))" <<'PYEOF'
import math, os, random, sys
refdir, out, k, seed = sys.argv[1:5]
random.seed(int(seed))
refs = [f for f in os.listdir(refdir) if f.endswith(".fa")]
pick = random.sample(refs, min(int(k), len(refs)))
abund = {f: math.e ** random.uniform(-2.3, 2.3) for f in pick}
with open(out, "w") as out_fh:
    for f in pick:
        with open(os.path.join(refdir, f)) as fh:
            for line in fh:
                if line.startswith(">"):
                    pid = f[:-3]
                    out_fh.write(">%s ab=%.4f\n" % (pid, abund[f]))
                else:
                    out_fh.write(line)
print("community:", len(pick), "plasmids")
PYEOF

      # Simulate Illumina reads (ART); fold = 150/250 etc. fragment size fixed
      art_illumina -ss HSXt -i "${comm}" \
        -l "${rl}" -f "${dep}" -m 350 -s 40 \
        -o reads_ -na
      # ART names reads_1.fq / reads_2.fq

      if [ "${ASSEMBLER}" = "megahit" ]; then
        megahit -1 reads_1.fq -2 reads_2.fq -o asm -t "${THREADS}" \
          --presets meta-sensitive
        ASM_FA="${rdir}/asm/final.contigs.fa"
      else
        metaspades.py -1 reads_1.fq -2 reads_2.fq -o asm -t "${THREADS}" -m 200
        ASM_FA="${rdir}/asm/contigs.fasta"
      fi

      # Map reads back and compute per-contig depth
      bowtie2-build -q "${ASM_FA}" asm_idx >/dev/null 2>&1
      bowtie2 -x asm_idx -1 reads_1.fq -2 reads_2.fq -p "${THREADS}" \
        --very-fast -S mapped.sam 2>/dev/null
      samtools view -bS mapped.sam | samtools sort -o mapped.bam
      samtools index mapped.bam
      jgi_summarize_bam_contig_depths --outputDepth depth.tsv mapped.bam

      # Predict and aggregate plasmids (geNomad `classify`; --clean-up is a
      # boolean flag and the positional order is INPUT OUTPUT_DIR DATABASE).
      GENOMAD_OUT="${rdir}/genomad_out"
      genomad classify --clean-up --threads "${THREADS}" \
        "${ASM_FA}" "${GENOMAD_OUT}" "${GENOMAD_DB:?set GENOMAD_DB}"
      plasrisk meta "${ASM_FA}" --genomad-bin "${GENOMAD_OUT}" \
        --coverage depth.tsv -o plasrisk_meta --min-contig 1000 \
        --complete-length 10000 --threads "${THREADS}"

      echo -e "${rep}\t${rl}\t${dep}\t${comm}\t${rdir}" >> "${GRID_LOG}"
      cd "${WORKDIR}"
    done
  done
done

# -----------------------------------------------------------------------------
# Match predicted bins to truth plasmids and score recovery.
# Matching rule: a bin matches a truth plasmid if >=50% of the truth plasmid's
# cargo gene names (ARG+BMG+VF; use column cargo_genes, else high_risk_genes)
# appear in the bin; choose the bin with max gene overlap. Plasmids with no
# cargo gene are reported as unmatchable, not as failures.
# -----------------------------------------------------------------------------
python3 - "${TRUTH_TABLE}" "${GRID_LOG}" <<'PYEOF'
import os, sys
import pandas as pd

truth_file, manifest_file = sys.argv[1], sys.argv[2]
truth = pd.read_csv(truth_file, sep="\t")
gcol = "cargo_genes" if "cargo_genes" in truth.columns else "high_risk_genes"
def geneset(x):
    return set(str(x).split(";")) - {"", "nan"}
truth["gs"] = truth[gcol].map(geneset)

rows = []
man = pd.read_csv(manifest_file, sep="\t")
for _, m in man.iterrows():
    bins = pd.read_csv(os.path.join(m["dir"], "plasrisk_meta",
                                   "meta_bin_scores.tsv"), sep="\t")
    def bin_genes(r):
        return (set(str(r["high_risk_genes"]).split(";"))
                | set(str(r.get("replicon", "")).split(";"))) - {""}
    bins["gs"] = bins.apply(bin_genes, axis=1)
    for _, t in truth.iterrows():
        if not t["gs"]:
            rows.append({
                "run": os.path.basename(m["dir"]),
                "plasmid_id": t["plasmid_id"], "matched": False,
                "S_truth": t["S_whole"], "grade_truth": t["grade"],
                "S_bin": None, "grade_bin": "", "confidence": "",
                "gene_recovery": None, "matchable": False})
            continue
        best, best_ov = None, 0
        for _, b in bins.iterrows():
            ov = len(t["gs"] & b["gs"])
            if ov > best_ov:
                best, best_ov = b, ov
        denom = max(1, len(t["gs"]))
        matched = best is not None and best_ov / denom >= 0.5
        rows.append({
            "run": os.path.basename(m["dir"]),
            "plasmid_id": t["plasmid_id"],
            "matched": matched,
            "S_truth": t["S_whole"],
            "grade_truth": t["grade"],
            "S_bin": best["S_lite"] if matched else None,
            "grade_bin": best["grade"] if matched else "",
            "confidence": best["confidence"] if matched else "",
            "gene_recovery": best_ov / denom,
            "matchable": True,
        })
res = pd.DataFrame(rows)
res.to_csv(os.path.join(os.path.dirname(manifest_file),
                        "fragment_grid_truth_matching.tsv"),
           sep="\t", index=False)
mm = res[res["matchable"]]
print(mm.groupby("run").agg(
    n_matchable=("plasmid_id", "count"),
    matched_pct=("matched", lambda x: 100 * x.mean()),
    mean_gene_recovery=("gene_recovery", "mean")))
PYEOF

echo "[grid] DONE. Summary: fragment_grid_truth_matching.tsv"
