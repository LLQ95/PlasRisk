#!/usr/bin/env python3
"""13_ablation_genomad_only.py

P0-3 ablation: geNomad-only versus geNomad + PlasRisk, run on the
existing prospective-cohort outputs (no re-assembly).

geNomad is the single recovery engine in both arms, so binary recovery
is identical by construction. PlasRisk adds risk ranking on top of the
recovered bins. The ablation quantifies that added value:
  - concordance of geNomad plasmid_score and PlasRisk S (Spearman);
  - ARG / high-risk-ARG enrichment in the top fraction selected by each
    ranking (top 1%, 5%, 10%);
  - enrichment at the PlasRisk grade A/B and grade >= C sets, compared
    with score-selected sets of the same size;
  - the eight whole-assembly critical-ARG contigs in the traveller
    isolate cohort (how the two arms classify each).

Usage:
  python3 13_ablation_genomad_only.py \
      --isolates results/traveler_isolates \
      --sewage   results/single_group \
      --cef      results/cefprozil \
      --out      results/ablation
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

TOP_FRACS = (0.01, 0.05, 0.10)

# canonical meta_bin_scores columns
CANON = ["bin_id", "plasmid_score", "n_ARG", "high_risk_genes",
         "S_lite", "S_frag4", "S_cargo3", "grade", "bin_origin"]


def spearman(a, b):
    try:
        from scipy.stats import spearmanr
        r = spearmanr(a, b)
        pv = float(r.pvalue)
        pstr = "<1e-300" if pv == 0 else ("%.3g" % pv)
        return round(float(r.statistic), 4), pstr
    except Exception:
        return "", ""


def _find(cols, exact, regexes):
    """Find a column: exact names first, then regex over lowered names."""
    for name in exact:
        if name in cols:
            return name
    lowered = {c.lower(): c for c in cols}
    for pat in regexes:
        rx = re.compile(pat, re.IGNORECASE)
        for low, orig in lowered.items():
            if rx.match(low):
                return orig
    return None


def _print_header(path, cols):
    sys.stderr.write("[schema] %s\n" % path)
    sys.stderr.write("         columns: %s\n" % " | ".join(cols))


def normalize_meta(path, d):
    """Map any historical meta_bin_scores schema onto canonical names."""
    cols = list(d.columns)
    mapping = {}
    arg = _find(cols,
                ["n_ARG", "n_ARG_genes", "n_arg", "n_args", "arg_count"],
                [r"^n_?arg(_genes|_count|_hits)?$", r"^arg_count$"])
    if arg is None:
        _print_header(path, cols)
        raise SystemExit(
            "[13] cannot locate an ARG-count column in %s\n"
            "    columns were: %s" % (path, " | ".join(cols)))
    mapping[arg] = "n_ARG"

    hr = _find(cols,
               ["high_risk_genes", "high_risk_args", "high_risk"],
               [r"^high[_-]?risk(_?(genes|args?))?$"])
    if hr is None:
        _print_header(path, cols)
        raise SystemExit(
            "[13] cannot locate a high-risk-genes column in %s\n"
            "    columns were: %s" % (path, " | ".join(cols)))
    mapping[hr] = "high_risk_genes"

    for src, dst in [("bin_id", "bin_id"), ("plasmid_score", "plasmid_score"),
                     ("S_lite", "S_lite"), ("S_frag4", "S_frag4"),
                     ("S_cargo3", "S_cargo3"), ("grade", "grade"),
                     ("bin_origin", "bin_origin")]:
        if src in cols:
            mapping[src] = dst

    d = d.rename(columns=mapping)
    # missing optional columns are created as empty
    for c in CANON:
        if c not in d.columns:
            d[c] = pd.NA
    return d[CANON + [c for c in d.columns if c not in CANON]], mapping


def read_isolates(root):
    scores_p = os.path.join(root, "benchmark_isolate_scores.tsv")
    if not os.path.exists(scores_p):
        return None, None
    sc = pd.read_csv(scores_p, sep="\t", low_memory=False)
    seq_col = next((c for c in ("seq_id", "seq_name", "contig")
                    if c in sc.columns), None)
    if seq_col is None:
        raise SystemExit("[13] isolate scores file has no seq_id column: %s"
                         % scores_p)
    rows = []
    for p in glob.glob(os.path.join(root, "genomad", "**",
                                    "*_plasmid_summary.tsv"),
                       recursive=True):
        try:
            d = pd.read_csv(p, sep="\t", low_memory=False)
        except Exception:
            continue
        if "seq_name" not in d.columns:
            continue
        d = d.rename(columns={"seq_name": seq_col})
        d = d[[seq_col] + [c for c in ("plasmid_score", "length", "topology")
                           if c in d.columns]]
        rows.append(d)
    g = (pd.concat(rows, ignore_index=True).drop_duplicates(seq_col)
         if rows else pd.DataFrame(columns=[seq_col]))
    return sc, g


def read_meta(root):
    frames = []
    seen_mappings = set()
    paths = sorted(glob.glob(os.path.join(root, "plasrisk", "*",
                                          "meta_bin_scores.tsv")))
    if not paths:
        # tolerate alternate layouts (e.g. extra cohort subdirectory)
        paths = sorted(glob.glob(os.path.join(root, "**",
                                              "meta_bin_scores.tsv"),
                                 recursive=True))
    for p in paths:
        try:
            d = pd.read_csv(p, sep="\t", low_memory=False)
        except Exception:
            continue
        d, mapping = normalize_meta(p, d)
        sid = os.path.basename(os.path.dirname(p))
        d["sample"] = sid
        key = tuple(sorted(mapping.items()))
        if key not in seen_mappings:
            seen_mappings.add(key)
            sys.stderr.write("[schema] %s -> %s\n" % (
                os.path.relpath(p, root),
                ", ".join("%s->%s" % (k, v) for k, v in sorted(mapping.items()))))
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    d = pd.concat(frames, ignore_index=True)
    if "bin_origin" in d.columns:
        d = d[d["bin_origin"].fillna("") == "predicted"].copy()
    return d


def rank_enrichment(df, score_col, n_arg, n_hr, tag):
    """ARG/high-risk capture at top fractions under one ranking."""
    d = df.dropna(subset=[score_col]).copy()
    out = []
    n = len(d)
    if n == 0:
        return out
    arg = d[n_arg].fillna(0).astype(float)
    hr = d[n_hr].fillna("").astype(str)
    for frac in TOP_FRACS:
        k = max(1, int(round(n * frac)))
        top = d.nlargest(k, score_col)
        top_arg = top[n_arg].fillna(0).astype(float)
        top_hr = top[n_hr].fillna("").astype(str)
        bg = float((arg > 0).mean())
        out.append({
            "ranking": tag, "top_frac": frac, "k": k,
            "arg_positive_in_set": int((top_arg > 0).sum()),
            "arg_frac_in_set": round(float((top_arg > 0).mean()), 4),
            "high_risk_in_set": int((top_hr != "").sum()),
            "arg_background_frac": round(bg, 4),
            "arg_enrichment": round(float((top_arg > 0).mean())
                                   / max(1e-9, bg), 2),
        })
    return out


def set_comparison(df, score_col, n_arg, n_hr, grade_col, tag):
    """Grade-set vs equal-size score-selected set."""
    out = []
    for glabel, mask in [("A_or_B", df[grade_col].isin(["A", "B"])),
                         ("C_or_above", df[grade_col].isin(["A", "B", "C"]))]:
        k = int(mask.sum())
        if k == 0:
            out.append({"comparison": tag, "grade_set": glabel,
                        "set_size": 0, "grade_arg": 0, "grade_high_risk": 0,
                        "score_arg": 0, "score_high_risk": 0})
            continue
        top = df.dropna(subset=[score_col]).nlargest(k, score_col)
        g_arg = df.loc[mask, n_arg].fillna(0).astype(float)
        g_hr = df.loc[mask, n_hr].fillna("").astype(str)
        out.append({
            "comparison": tag, "grade_set": glabel, "set_size": k,
            "grade_arg": int((g_arg > 0).sum()),
            "grade_high_risk": int((g_hr != "").sum()),
            "score_arg": int((top[n_arg].fillna(0).astype(float) > 0).sum()),
            "score_high_risk": int((top[n_hr].fillna("").astype(str) != "").sum()),
        })
    return out


def analyze_meta(df, tag):
    if df.empty:
        return [], []
    s = df["S_lite"].combine_first(df["S_frag4"]).combine_first(df["S_cargo3"])
    df = df.assign(PlasRisk_S=s)
    rows = []
    pair = df.dropna(subset=["plasmid_score", "PlasRisk_S"])
    rho, p = spearman(pair["plasmid_score"], pair["PlasRisk_S"]) if len(pair) else ("", "")
    rows.append({"cohort": tag, "n_predicted_bins": len(df),
                 "spearman_score_vs_S": rho, "spearman_p": p})
    rank_rows, set_rows = [], []
    use = df.dropna(subset=["plasmid_score"])
    for r in rank_enrichment(use, "plasmid_score", "n_ARG",
                             "high_risk_genes", "geNomad"):
        r["cohort"] = tag
        rank_rows.append(r)
    use2 = df.dropna(subset=["PlasRisk_S"])
    for r in rank_enrichment(use2, "PlasRisk_S", "n_ARG",
                             "high_risk_genes", "PlasRisk"):
        r["cohort"] = tag
        rank_rows.append(r)
    for r in set_comparison(df, "plasmid_score", "n_ARG", "high_risk_genes",
                            "grade", tag):
        set_rows.append(r)
    return rank_rows, set_rows + [{"summary": rows}]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isolates", default="")
    ap.add_argument("--sewage", default="")
    ap.add_argument("--cef", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    all_rank, all_set, overview = [], [], []

    # metagenome cohorts
    for tag, root in [("sewage_110", args.sewage), ("cefprozil_72", args.cef)]:
        if root and os.path.isdir(root):
            df = read_meta(root)
            if df.empty:
                sys.stderr.write(
                    "[13] WARNING: no meta_bin_scores.tsv found under %s; "
                    "cohort %s skipped\n" % (root, tag))
            if not df.empty:
                r, s = analyze_meta(df, tag)
                all_rank.extend(r)
                ov = s.pop()
                overview.extend(ov["summary"])
                all_set.extend(s)
                df.to_csv(os.path.join(args.out,
                          "ablation_bins_%s.tsv" % tag), sep="\t", index=False)

    # isolate cohort
    if args.isolates and os.path.isdir(args.isolates):
        sc, g = read_isolates(args.isolates)
        if sc is not None and not g.empty:
            df = sc.merge(g, on="seq_id", how="left")
            df["high_risk"] = df["high_risk_genes"].fillna("").astype(str)
            df["ARG"] = df["n_ARG"].fillna(0).astype(float)
            pair = df.dropna(subset=["plasmid_score", "S_norm"])
            rho, p = spearman(pair["plasmid_score"], pair["S_norm"]) if len(pair) else ("", "")
            overview.append({"cohort": "isolates_189",
                             "n_plasmid_contigs": len(df),
                             "spearman_score_vs_S": rho, "spearman_p": p})
            for frac in TOP_FRACS:
                k = max(1, int(round(len(df) * frac)))
                gn = df.nlargest(k, "plasmid_score")
                pr = df.nlargest(k, "S_norm")
                bg = float((df["ARG"] > 0).mean())
                for sub, rname in [(gn, "geNomad"), (pr, "PlasRisk")]:
                    all_rank.append({
                        "cohort": "isolates_189", "ranking": rname,
                        "top_frac": frac, "k": k,
                        "arg_positive_in_set": int((sub["ARG"] > 0).sum()),
                        "arg_frac_in_set": round(float((sub["ARG"] > 0).mean()), 4),
                        "high_risk_in_set": int((sub["high_risk"] != "").sum()),
                        "arg_background_frac": round(bg, 4),
                        "arg_enrichment": round(float((sub["ARG"] > 0).mean())
                                                / max(1e-9, bg), 2)})
            for glabel, mask in [("A_or_B", df["grade"].isin(["A", "B"])),
                                 ("C_or_above", df["grade"].isin(["A", "B", "C"]))]:
                k = int(mask.sum())
                if k:
                    top = df.nlargest(k, "plasmid_score")
                    g_arg = df.loc[mask, "ARG"]
                    g_hr = df.loc[mask, "high_risk"]
                    all_set.append({
                        "comparison": "isolates_189", "grade_set": glabel,
                        "set_size": k,
                        "grade_arg": int((g_arg > 0).sum()),
                        "grade_high_risk": int((g_hr != "").sum()),
                        "score_arg": int((top["ARG"] > 0).sum()),
                        "score_high_risk": int((top["high_risk"] != "").sum())})
            # the 8 whole-assembly critical-ARG contigs
            crit_p = os.path.join(args.isolates, "critical_gene_hits.tsv")
            if os.path.exists(crit_p):
                crit = pd.read_csv(crit_p, sep="\t", dtype=str, low_memory=False)
                ccol = next((c for c in ("seq_id", "contig", "contig_id")
                             if c in crit.columns), None)
                if ccol:
                    m = df.set_index("seq_id")
                    crit["geNomad_plasmid_call"] = crit[ccol].map(
                        lambda x: "called" if x in m.index else "NOT_called")
                    crit["geNomad_plasmid_score"] = crit[ccol].map(
                        lambda x: m.loc[x, "plasmid_score"] if x in m.index else "")
                    crit["PlasRisk_grade"] = crit[ccol].map(
                        lambda x: m.loc[x, "grade"] if x in m.index else "")
                    crit.to_csv(os.path.join(args.out,
                               "ablation_critical_contigs.tsv"),
                               sep="\t", index=False)

    pd.DataFrame(overview).to_csv(os.path.join(args.out,
        "ablation_overview.tsv"), sep="\t", index=False)
    pd.DataFrame(all_rank).to_csv(os.path.join(args.out,
        "ablation_top_fraction_enrichment.tsv"), sep="\t", index=False)
    pd.DataFrame(all_set).to_csv(os.path.join(args.out,
        "ablation_grade_vs_score_sets.tsv"), sep="\t", index=False)

    print("=== ablation overview ===")
    print(pd.DataFrame(overview).to_string(index=False))
    print("\n=== top-fraction enrichment ===")
    print(pd.DataFrame(all_rank).to_string(index=False))
    print("\n=== grade set vs equal-size score set ===")
    print(pd.DataFrame(all_set).to_string(index=False))
    print("\nwritten to", args.out)


if __name__ == "__main__":
    sys.exit(main())
