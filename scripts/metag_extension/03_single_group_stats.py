#!/usr/bin/env python3
"""
03_single_group_stats.py
Single-habitat-group descriptive statistics for the preassembled benchmark
(used for the 110 global wastewater metagenomes, Hendriksen et al. 2019,
when the 24 plasmidome samples are not included).

Cross-group trend tests (Spearman/Kruskal/Cliff's delta) are undefined for a
single group and produce NaN; this script does not compute them. It reports:

  - per-metric descriptive statistics with a 95% CI of the mean
  - sample-level prevalence (fraction of samples with a positive value)
    with Wilson 95% CI for each percentage metric
  - country-level summaries (country parsed from the accession prefix,
    e.g. CAN-22_1 -> CAN; sites counted after stripping _<rep> suffix)
  - Spearman correlation matrix among metrics
  - optional independent-truth correlation (--external-truth)

Usage:
  python3 03_single_group_stats.py benchmark_sample_summary.tsv outdir/
      [--group wastewater]
      [--external-truth truth.tsv]

external-truth TSV (header required, tab-separated), one of:
  country <TAB> external_metric      (per-country median comparison)
  accession <TAB> external_metric    (per-sample comparison)
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

METRIC_NAMES = [
    "CPRI_mean", "CPRI_p90", "pct_grade_AB",
    "critical_ARG_load_per_Mb", "pct_ARG_positive",
    "pct_high_risk_ARG",
]
# Predicted-plasmid-bin metrics (plasrisk >= 1.7.3): caller-predicted bins
# only, free of whole-community chromosome-fragment dilution.
PRED_METRIC_NAMES = [
    "pred_CPRI_mean", "pred_CPRI_p90", "pred_pct_grade_AB",
    "pred_critical_ARG_load_per_Mb", "pred_pct_ARG_positive",
    "pred_pct_high_risk_ARG",
]
PRIMARY_METRICS = METRIC_NAMES + PRED_METRIC_NAMES
# Percentage metrics for which a sample-level positive prevalence is meaningful.
PERCENT_METRICS = [
    "pct_grade_AB", "pct_ARG_positive", "pct_high_risk_ARG",
    "pred_pct_grade_AB", "pred_pct_ARG_positive",
    "pred_pct_high_risk_ARG",
]
# Bin-count columns carried through to the per-sample annotated table.
COUNT_COLS = ["n_plasmid_bins", "n_predicted_bins", "n_unbinned_bins",
              "n_scored_bins", "pred_n_scored_bins", "n_graded_bins",
              "pred_n_graded_bins"]


def wilson_ci(k, n, z=1.96):
    """Wilson score interval for a binomial proportion k/n."""
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / den
    return p, max(0.0, centre - half), min(1.0, centre + half)


def country_of(acc):
    m = re.match(r"^([A-Za-z]{3})", str(acc))
    return m.group(1).upper() if m else "UNK"


def site_of(acc):
    # CAN-22_1 / CAN-22a_3 -> CAN-22 / CAN-22a ; KHM-21 stays KHM-21.
    return re.sub(r"_\d+$", "", str(acc))


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


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Single-group descriptive statistics (no gradient tests).")
    ap.add_argument("infile", help="benchmark_sample_summary.tsv")
    ap.add_argument("outdir", help="directory for output tables")
    ap.add_argument(
        "--group", default=None,
        help="habitat_group to analyze; required if the table contains "
             "more than one group (default: the sole group present)")
    ap.add_argument(
        "--external-truth", default=None,
        help="optional TSV with independent truth (country,external_metric "
             "or accession,external_metric)")
    args = ap.parse_args(argv)

    if not os.path.exists(args.infile):
        sys.exit("[error] input not found: %s" % args.infile)
    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.infile, sep="\t")
    if "habitat_group" not in df.columns:
        sys.exit("[error] input has no habitat_group column: %s" % args.infile)

    groups = sorted(df["habitat_group"].dropna().astype(str).unique().tolist())
    if args.group:
        df = df[df["habitat_group"] == args.group].copy()
        if df.empty:
            sys.exit("[error] no rows for group '%s'; present: %s"
                     % (args.group, groups))
    elif len(groups) == 1:
        df = df[df["habitat_group"] == groups[0]].copy()
    else:
        sys.exit("[error] table contains %d groups %s; rerun with --group "
                 "to select one." % (len(groups), groups))
    group = df["habitat_group"].iloc[0]
    df["country"] = df["accession"].map(country_of)
    df["site"] = df["accession"].map(site_of)

    metrics = [m for m in PRIMARY_METRICS if m in df.columns]
    if not metrics:
        sys.exit("[error] none of the expected metric columns are present: %s"
                 % PRIMARY_METRICS)
    print("group: %s | samples: %d | countries: %d | sites: %d"
          % (group, len(df), df["country"].nunique(), df["site"].nunique()))
    if "n_predicted_bins" in df.columns:
        print("median bins/sample: predicted=%.0f unbinned=%.0f"
              % (df["n_predicted_bins"].median(),
                 df["n_unbinned_bins"].median()
                 if "n_unbinned_bins" in df.columns else float("nan")))

    # 1) Per-metric descriptive statistics with 95% CI of the mean.
    desc_rows = []
    for m in metrics:
        x = pd.to_numeric(df[m], errors="coerce").dropna()
        n = len(x)
        mean = x.mean() if n else np.nan
        sd = x.std(ddof=1) if n > 1 else np.nan
        if n > 1:
            tcrit = stats.t.ppf(0.975, n - 1)
            half = tcrit * sd / np.sqrt(n)
            lo, hi = mean - half, mean + half
        else:
            lo = hi = np.nan
        desc_rows.append({
            "metric": m, "n": n, "mean": mean, "sd": sd,
            "ci95_low": lo, "ci95_high": hi,
            "median": x.median() if n else np.nan,
            "q1": x.quantile(0.25) if n else np.nan,
            "q3": x.quantile(0.75) if n else np.nan,
            "min": x.min() if n else np.nan,
            "max": x.max() if n else np.nan,
        })
    desc = pd.DataFrame(desc_rows)
    desc.to_csv(os.path.join(args.outdir, "single_group_descriptive.tsv"),
                sep="\t", index=False)

    # 2) Sample-level prevalence of positive values (Wilson 95% CI).
    prev_rows = []
    n_total = len(df)
    for m in PERCENT_METRICS:
        if m not in df.columns:
            continue
        x = pd.to_numeric(df[m], errors="coerce")
        k = int((x > 0).sum())
        n = int(x.notna().sum())
        p, lo, hi = wilson_ci(k, n)
        prev_rows.append({
            "metric": m, "samples_positive": k, "samples_with_data": n,
            "prevalence": p, "wilson_ci95_low": lo, "wilson_ci95_high": hi,
        })
    prev = pd.DataFrame(prev_rows)
    prev.to_csv(os.path.join(args.outdir, "sample_prevalence_wilson.tsv"),
                sep="\t", index=False)

    # 3) Country-level summaries (replicate samples/site collapsed to median).
    agg = {m: "median" for m in metrics}
    agg["site"] = pd.Series.nunique
    by_country = (df.groupby("country")
                    .agg({**agg, "accession": "count"})
                    .rename(columns={"accession": "n_samples"}))
    by_country = by_country.sort_values("n_samples", ascending=False)
    by_country.to_csv(os.path.join(args.outdir, "by_country.tsv"),
                      sep="\t")

    # 4) Spearman correlation matrix among metrics.
    corr = df[metrics].apply(pd.to_numeric, errors="coerce").corr(
        method="spearman")
    corr.to_csv(os.path.join(args.outdir, "metric_correlation.tsv"),
                sep="\t")

    # 5) Per-sample annotated table (country, site, bin counts).
    annot_cols = ["accession", "habitat_group", "country", "site"]
    annot_cols += [c for c in COUNT_COLS if c in df.columns]
    annot_cols += metrics
    df[annot_cols].to_csv(
        os.path.join(args.outdir, "samples_annotated.tsv"),
        sep="\t", index=False)

    # 6) Optional independent-truth validation.
    if args.external_truth:
        truth = pd.read_csv(args.external_truth, sep="\t")
        if "accession" in truth.columns:
            merged = df.merge(truth, on="accession", how="inner")
        elif "country" in truth.columns:
            med = df.groupby("country")[metrics].median()
            merged = med.merge(
                truth.set_index("country")["external_metric"],
                left_index=True, right_index=True, how="inner")
        else:
            sys.exit("[error] truth file needs country,external_metric or "
                     "accession,external_metric columns")
        if merged.empty:
            sys.exit("[error] truth file matched no samples/countries")
        rows = []
        for m in metrics:
            d = merged[[m, "external_metric"]].dropna()
            if len(d) < 4:
                continue
            rho, p = stats.spearmanr(d[m], d["external_metric"])
            rows.append({"metric": m, "n": len(d),
                         "spearman_rho_vs_truth": rho, "p": p})
        ev = pd.DataFrame(rows)
        if not ev.empty:
            ev["p_BH"] = bh_adjust(ev["p"].fillna(1))
        ev.to_csv(os.path.join(args.outdir, "external_validation.tsv"),
                  sep="\t", index=False)
        print("\n=== independent truth correlation (n=%d matched units) ==="
              % len(merged))
        print(ev.to_string(index=False) if not ev.empty else "(n<4)")

    print("\n=== descriptive statistics ===")
    print(desc.to_string(index=False))
    print("\n=== sample-level prevalence (Wilson 95% CI) ===")
    print(prev.to_string(index=False) if not prev.empty else "(none)")
    print("\nwrote", args.outdir)


if __name__ == "__main__":
    main()
