#!/usr/bin/env python3
"""
03_gradient_stats.py
Tier A statistics for the anthropogenic-intensity habitat gradient.

Usage:
  python3 03_gradient_stats.py benchmark_sample_summary.tsv outdir/

The input is the table produced by 02_public_metag_benchmark.sh. Habitat
groups are mapped to a prespecified anthropogenic-intensity order (edit the
ORDER below; it must be fixed before looking at results).

Tests:
  - Spearman trend of each metric against the ordered group index
  - Kruskal-Wallis across habitat groups
  - paired Wilcoxon for samples sharing a pair_id
  - BH adjustment across metrics
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

# Prespecified order, low -> high anthropogenic antibiotic pressure.
# Edit to match the habitat labels used in accessions.tsv BEFORE running.
ORDER = [
    "healthy_human_gut",
    "receiving_water_soil",
    "livestock_manure",
    "antibiotic_exposed_human_gut",
    "hospital_influent",
]
METRICS = [
    "CPRI_mean", "CPRI_p90", "pct_grade_AB",
    "critical_ARG_load_per_Mb", "pct_ARG_positive",
    "pct_high_risk_ARG",
]


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def cliffs_delta(x, gidx):
    """Effect size: probability a sample in a higher group has a higher
    value minus the reverse (rank-biserial analogue for the trend)."""
    x = np.asarray(x)
    g = np.asarray(gidx)
    num, den = 0.0, 0.0
    for i in range(len(x)):
        for j in range(len(x)):
            if g[j] > g[i]:
                den += 1
                num += np.sign(x[j] - x[i])
    return num / den if den else np.nan


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Tier A anthropogenic-intensity gradient statistics.")
    ap.add_argument(
        "infile",
        help="benchmark_sample_summary.tsv produced by 02_public_metag_benchmark.sh")
    ap.add_argument("outdir", help="directory for output tables")
    ap.add_argument(
        "--order", default=",".join(ORDER),
        help="comma-separated habitat_group labels, low->high intensity; "
             "must match the labels in accessions.tsv (default: prespecified "
             "5-level order)")
    args = ap.parse_args(argv)
    order = [g.strip() for g in args.order.split(",") if g.strip()]

    if not os.path.exists(args.infile):
        sys.exit(
            "[error] input not found: %s\n"
            "        03 reads benchmark_sample_summary.tsv, which is written "
            "to the workdir only after 02_public_metag_benchmark.sh has run its "
            "process stage (fastp->megahit->genomad->plasrisk meta). Finish that "
            "first, then rerun 03." % args.infile)
    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.infile, sep="\t")
    if "habitat_group" not in df.columns:
        sys.exit("[error] input has no habitat_group column: %s" % args.infile)
    present = sorted(df["habitat_group"].dropna().astype(str).unique().tolist())
    df = df[df["habitat_group"].isin(order)].copy()
    if df.empty:
        sys.exit(
            "[error] no rows match the habitat order.\n"
            "        order given : %s\n"
            "        present     : %s\n"
            "        rerun with --order matching the present labels, e.g. "
            "--order wastewater,livestock" % (order, present))
    df["group_idx"] = df["habitat_group"].map(
        {g: i for i, g in enumerate(order)})
    print("n samples by group:")
    print(df.groupby("habitat_group").size().reindex(order))

    trend_rows = []
    for metric in METRICS:
        if metric not in df.columns:
            continue
        d = df.dropna(subset=[metric])
        if len(d) < 4:
            continue
        rho, p_spear = stats.spearmanr(d["group_idx"], d[metric])
        groups = [g[metric].to_numpy() for _, g in d.groupby("group_idx")
                  if g[metric].notna().any()]
        if len(groups) >= 2:
            H, p_kw = stats.kruskal(*groups)
        else:
            H, p_kw = np.nan, np.nan
        trend_rows.append({
            "metric": metric, "n": len(d),
            "spearman_rho": rho, "spearman_p": p_spear,
            "kruskal_H": H, "kruskal_p": p_kw,
            "cliffs_delta_vs_order": cliffs_delta(d[metric], d["group_idx"]),
        })
    trend = pd.DataFrame(trend_rows)
    if not trend.empty:
        trend["spearman_p_BH"] = bh_adjust(trend["spearman_p"])
        trend["kruskal_p_BH"] = bh_adjust(trend["kruskal_p"].fillna(1))
    trend.to_csv(os.path.join(args.outdir, "gradient_trend_tests.tsv"),
                 sep="\t", index=False)

    # Paired comparisons (pair_id must be shared by exactly two samples)
    pair_rows = []
    if "pair_id" in df.columns:
        paired = df[df["pair_id"].astype(str) != "-"]
        for pid, pp in paired.groupby("pair_id"):
            if len(pp) != 2:
                continue
            for metric in METRICS:
                if metric not in pp.columns or pp[metric].isna().any():
                    continue
                a, b = pp[metric].to_numpy()
                try:
                    stat, p = stats.wilcoxon(a, b)
                except ValueError:
                    stat, p = np.nan, np.nan
                pair_rows.append({
                    "pair_id": pid, "metric": metric,
                    "value_first": a, "value_second": b,
                    "delta": b - a, "wilcoxon_p": p,
                })
    pairs = pd.DataFrame(pair_rows)
    if not pairs.empty:
        pairs["wilcoxon_p_BH"] = bh_adjust(pairs["wilcoxon_p"].fillna(1))
    pairs.to_csv(os.path.join(args.outdir, "paired_comparisons.tsv"),
                 sep="\t", index=False)

    # Group medians/IQRs for figure construction
    desc = df.groupby("habitat_group")[METRICS].agg(
        ["median", lambda q: q.quantile(0.25), lambda q: q.quantile(0.75),
         "count"])
    desc.to_csv(os.path.join(args.outdir, "group_descriptive.tsv"), sep="\t")

    print("\n=== trend tests ===")
    print(trend.to_string(index=False) if not trend.empty else "(none)")
    print("\nwrote", args.outdir)


if __name__ == "__main__":
    main()
