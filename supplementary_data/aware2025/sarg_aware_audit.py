# -*- coding: utf-8 -*-
"""Recompute S_ARG under (a) old regex AWaRe mapping and (b) 2025 AWaRe mapping;
compare with locked component-matrix S_ARG; quantify discrimination impact."""
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score

MASTER = r"C:\Users\Administrator\Doubao\chats\2026-08-15\new-chat\pipdb_pipeline\results_run\psc_master.tsv"
LOCKED_CM = r"C:\Users\Administrator\Downloads\data\salmonella\plasmid_risk_assessment\revision_2026.9.21\data_release\plasrisk_component_matrix_792964x10.csv"
SPLIT = r"C:\Users\Administrator\Doubao\chats\2026-08-15\new-chat\pipdb_pipeline\results_run\tables\locked_split_ids.csv.gz"
OUTDIR = r"C:\Users\Administrator\Doubao\chats\2026-08-30\new-chat"

# ---------- OLD classifier (exact replica of pipdb_02 AWARE_RULES) ----------
AWARE_RULES = [
    (r"penem|carbapenem|monobactam", "Watch", 1),
    (r"cephalosporin", "Watch", 0),
    (r"penam|beta-lactam|β-lactam", "Access", 0),
    (r"aminoglycoside", "Watch", 0),
    (r"macrolide", "Watch", 0),
    (r"tetracycline", "Access", 0),
    (r"glycylcycline", "Watch", 1),
    (r"phenicol|chloramphenicol", "Access", 0),
    (r"sulfonamide|sulfone", "Access", 0),
    (r"diaminopyrimidine|trimethoprim", "Access", 0),
    (r"fluoroquinolone|quinolone", "Watch", 0),
    (r"glycopeptide", "Watch", 1),
    (r"oxazolidinone|linezolid", "Reserve", 1),
    (r"polymyxin|colistin", "Reserve", 1),
    (r"fosfomycin|phosphonic", "Reserve", 1),
    (r"rifamycin|rifampin", "Watch", 0),
    (r"nitroimidazole|nitrofuran", "Access", 0),
    (r"lincosamide|streptogramin", "Access", 0),
    (r"disinfecting|antiseptic|biocide", "Watch", 0),
    (r"nucleoside|peptide|lipopeptide", "Access", 0),
]
AWARE_W = {"Reserve": 3.0, "Watch": 2.0, "Access": 1.0, "unclassified": 1.0}

# ---------- NEW 2025 mapping: per-token category, then take the highest ----------
# category numeric: Access=1 Watch=2 Reserve=3 ; tokens for non-antibiotics excluded
TOKEN_CAT = [
    (r"carbapenem", 2), (r"\bpenem\b", 3), (r"monobactam", 3),
    (r"glycylcycline", 3),
    (r"cephalosporin|cephamycin", 2),
    (r"penam|beta-lactam|β-lactam", 1),
    (r"aminoglycoside", 2), (r"macrolide", 2),
    (r"tetracycline", 2),
    (r"phenicol|chloramphenicol", 1),
    (r"sulfonamide|sulfone", 1),
    (r"diaminopyrimidine|trimethoprim", 1),
    (r"fluoroquinolone|quinolone", 2),
    (r"glycopeptide", 2),
    (r"oxazolidinone|linezolid", 3),
    (r"polymyxin|colistin", 3),
    (r"fosfomycin|phosphonic", 3),
    (r"rifamycin|rifampin", 2),
    (r"nitroimidazole|nitrofuran", 1),
    (r"lincosamide", 2), (r"streptogramin", 3),
    (r"lipopeptide", 3),
    (r"nucleoside|peptide", 1),
]
CATNAME = {1: "Access", 2: "Watch", 3: "Reserve"}

HR_PATTERNS = [
    r"KPC|NDM|VIM|IMP|OXA-?48|OXA-?23|OXA-?24|OXA-?58|GES|SME|IMI|SPM|SIM|GIM",
    r"mcr-?\d", r"van[A-Z]\b", r"optra|^cfr|poxta", r"tet\(?X\)?|tetX",
    r"qnr[A-Z]?|aac\(6'\)-Ib-cr|qepa|oqxa", r"arma|rmt[A-H]|npma",
    r"fos[A-Z]?\d*", r"mec[A-C]|pvl",
]

def hr_genes(aro):
    if not isinstance(aro, str) or aro in ("\\N", ""):
        return [], 0
    genes = [g.strip() for g in aro.split(",") if g.strip()]
    n_hr = 0
    for g in genes:
        for hp in HR_PATTERNS:
            if re.search(hp, g, re.IGNORECASE):
                n_hr += 1
                break
    return genes, n_hr

def old_cat(drugclass):
    classes = str(drugclass).lower() if isinstance(drugclass, str) else ""
    aware, last = "unclassified", 0
    for pat, aw, lr in AWARE_RULES:
        if re.search(pat, classes):
            aware, last = aw, lr
            break
    return AWARE_W[aware] * (1.5 if last else 1.0), aware + ("*" if last else "")

def new_cat(drugclass):
    classes = str(drugclass).lower() if isinstance(drugclass, str) else ""
    tokens = [t.strip() for t in classes.split(",") if t.strip()]
    best = 1  # unclassified default Access weight, as in old code
    hit = False
    for tok in tokens:
        for pat, cat in TOKEN_CAT:
            if re.search(pat, tok):
                best = max(best, cat)
                hit = True
                break
    return float(best), CATNAME[best] + ("" if hit else "(unclassified)")

print("loading master ...")
df = pd.read_csv(MASTER, sep="\t", low_memory=False)
print("loading locked component matrix ...")
cm = pd.read_csv(LOCKED_CM)
print(cm.columns.tolist(), cm.shape)

res = df["aro_name"].apply(hr_genes)
df["n_hr"] = [r[1] for r in res]
old = df["drugclass"].apply(old_cat)
new = df["drugclass"].apply(new_cat)
df["wcat_old"] = [r[1] for r in old]
df["wcat_new"] = [r[1] for r in new]

def build_weight(base):
    w = base.copy()
    aro = df["aro_name"].fillna("\\N")
    for i, a in enumerate(aro.values):
        if not isinstance(a, str) or a in ("\\N", ""):
            w[i] = 0.0
            continue
        for g in [x.strip() for x in a.split(",") if x.strip()]:
            for hp in HR_PATTERNS:
                if re.search(hp, g, re.IGNORECASE):
                    w[i] *= 1.3
                    break
    return w

w_old = build_weight(np.array([r[0] for r in old], dtype=float))
w_new = build_weight(np.array([r[0] for r in new], dtype=float))

def sarg(w):
    cap = np.quantile(w, 0.99)
    return np.log10(1 + np.clip(w, None, cap)) / np.log10(1 + cap)

df["S_ARG_old"] = sarg(w_old)
df["S_ARG_new"] = sarg(w_new)
df["S_ARG_locked"] = cm["S_ARG"].values

y = (df["n_hr"] > 0).astype(int)
print("\nhigh-risk ARG positives:", int(y.sum()), f"({y.mean()*100:.2f}%)")

print("\n=== category assignment cross-tab (PSC counts) ===")
ct = pd.crosstab(df["wcat_old"], df["wcat_new"])
print(ct.to_string())

def stats(a, b, label):
    d = np.abs(a - b)
    print(f"\n--- {label} ---")
    print(f"  Spearman rho = {spearmanr(a, b).statistic:.6f}")
    print(f"  Pearson  r   = {pearsonr(a, b)[0]:.6f}")
    print(f"  mean |dS|    = {d.mean():.6f} ; median = {np.median(d):.6f} ; max = {d.max():.4f}")
    for t in [0.001, 0.01, 0.05, 0.1]:
        print(f"  PSCs with |dS| >= {t}: {(d >= t).sum():,} ({(d >= t).mean()*100:.2f}%)")

stats(df["S_ARG_locked"], df["S_ARG_old"], "locked vs OLD replica")
stats(df["S_ARG_locked"], df["S_ARG_new"], "locked vs NEW 2025")
stats(df["S_ARG_old"], df["S_ARG_new"], "OLD vs NEW")

print("\n=== AUC for high-risk ARG endpoint (all 792,964) ===")
for col in ["S_ARG_old", "S_ARG_new", "S_ARG_locked"]:
    print(f"  {col:14s} AUC = {roc_auc_score(y, df[col]):.6f}")

# locked test subset
sp = pd.read_csv(SPLIT)
te = sp.loc[sp["set"] == "locked_test", "id"].values - 1
yt = y.iloc[te]
print(f"\n=== AUC on locked test (n={len(te):,}, pos={yt.sum()}) ===")
for col in ["S_ARG_old", "S_ARG_new", "S_ARG_locked"]:
    print(f"  {col:14s} AUC = {roc_auc_score(yt, df[col].iloc[te]):.6f}")

print("\nmeans: old %.4f new %.4f locked %.4f" % (
    df["S_ARG_old"].mean(), df["S_ARG_new"].mean(), df["S_ARG_locked"].mean()))

ct.to_csv(OUTDIR + r"\tab_aware_mapping_crosstab.csv")
df[["id", "S_ARG_old", "S_ARG_new", "S_ARG_locked", "wcat_old", "wcat_new"]].to_csv(
    OUTDIR + r"\tab_sarg_aware_audit.csv", index=False)
print("\nsaved.")
