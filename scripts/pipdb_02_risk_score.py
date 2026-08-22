#!/usr/bin/env python3
"""
pipdb_02_risk_score.py
======================
PIPdb-adapted composite risk score + Q1-Q4 grading (8-dimension prototype).
The final 10-dimension PlasRisk weights are in pipdb_15/16/17 R scripts.

Usage:
  python pipdb_02_risk_score.py --config ../config.yaml
"""
import argparse, os, sys, re
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def load_config(p):
    import yaml
    with open(p) as fh: cfg = yaml.safe_load(fh)
    root = os.path.dirname(os.path.abspath(p))
    cfg["_root"]=root; cfg["_db"]=os.path.join(root,cfg["paths"]["pipdb_dir"])
    cfg["_out"]=os.path.join(root,cfg["paths"]["out_dir"])
    cfg["_ref"]=os.path.join(root,cfg["paths"]["ref_dir"])
    os.makedirs(cfg["_out"], exist_ok=True)
    return cfg

def log(m): print(f"[{pd.Timestamp.now():%H:%M:%S}] {m}", flush=True)

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
AWARE_W = {"Reserve":3.0, "Watch":2.0, "Access":1.0, "unclassified":1.0}

HR_PATTERNS = [
    r"KPC|NDM|VIM|IMP|OXA-?48|OXA-?23|OXA-?24|OXA-?58|GES|SME|IMI|SPM|SIM|GIM",
    r"mcr-?\d", r"van[A-Z]\b", r"optra|^cfr|poxta", r"tet\(?X\)?|tetX",
    r"qnr[A-Z]?|aac\(6'\)-Ib-cr|qepa|oqxa", r"arma|rmt[A-H]|npma",
    r"fos[A-Z]?\d*", r"mec[A-C]|pvl",
]

def arg_hazard(aro, drugclass, who, n_who):
    if not isinstance(aro,str) or aro in ("\\N",""):
        return 0.0, 0, 0
    genes = [g.strip() for g in aro.split(",") if g.strip()]
    classes = str(drugclass).lower() if isinstance(drugclass,str) else ""
    aware = "unclassified"; last = 0
    for pat, aw, lr in AWARE_RULES:
        if re.search(pat, classes):
            aware = aw; last = lr; break
    w = AWARE_W[aware] * (1.5 if last else 1.0)
    n_hr = 0
    for g in genes:
        gl = g.lower()
        for hp in HR_PATTERNS:
            if re.search(hp, gl, re.IGNORECASE):
                n_hr += 1; w *= 1.3; break
    if isinstance(who,str) and who not in ("\\N",""):
        w *= 1.2
    return w, len(genes), n_hr

def replicon_prior(rep):
    if pd.isna(rep): return 0.30
    r = str(rep)
    high = ["IncFII","IncFIA","IncFIB","IncFIC","IncA/C","IncA/C2","IncN","IncHI1",
            "IncHI2","IncX1","IncX3","IncX4","IncL/M","IncI1","IncI2","IncR","IncP",
            "IncP1","IncW","IncY","IncH","ColRNAI","ColKP3"]
    med  = ["Col156","Col440I","Col440II","Col(MG828)","Col(BS512)","IncQ1","IncQ2",
            "IncU","IncT","IncB/O/K/Z"]
    if any(r.startswith(h) or h in r for h in high): return 0.90
    if any(r.startswith(m) or m in r for m in med): return 0.50
    return 0.30

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="../config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)
    W = cfg["weights"]; G = cfg["grading"]

    log("loading psc_master ...")
    df = pd.read_csv(os.path.join(cfg["_out"],"psc_master.tsv"), sep="\t", low_memory=False)
    log(f"  {len(df):,} PSCs")

    log("computing ARG hazard ...")
    res = df.apply(lambda r: arg_hazard(r["aro_name"], r["drugclass"],
                                        r.get("arg_WHO","\\N"), r.get("n_who_arg",0)), axis=1)
    df["arg_weight"] = [r[0] for r in res]
    df["n_arg_computed"] = [r[1] for r in res]
    df["n_high_risk_arg"] = [r[2] for r in res]
    cap = df["arg_weight"].quantile(0.99)
    df["S_ARG"] = (np.log10(1+df["arg_weight"].clip(upper=cap))/np.log10(1+cap)).fillna(0)

    mob_map = {"conjugative_complete":1.0, "conjugative_likely":0.85,
               "mobilizable":0.6, "non-mobilizable":0.1}
    df["mob_base"] = df["mobility_class"].map(mob_map)
    df["mob_bonus"] = np.where(df["has_integron"], 0.1, 0.0)
    is_cap = 2.0
    df["S_MOB"] = (0.7*(df["mob_base"]+df["mob_bonus"]).clip(upper=1.0) +
                   0.3*(df["is_density_per_kb"].clip(upper=is_cap)/is_cap)).clip(0,1)

    hc_map = {"escape":1.0, "who_critical":0.9, "pathogen_other":0.6, "unknown":0.4}
    df["host_base"] = df["host_class"].map(hc_map)
    df["host_human_bonus"] = np.where(df["host_human"]|df["habitat_human"], 0.05, 0.0)
    vf_cap = 10
    df["S_HOST"] = (0.8*(df["host_base"]+df["host_human_bonus"]).clip(upper=1.0) +
                    0.2*(df["n_vf"].clip(upper=vf_cap)/vf_cap)).clip(0,1)

    df["S_REP"] = df["replicon_primary"].map(replicon_prior)

    size = pd.to_numeric(df["length_avg"], errors="coerce").fillna(1000).clip(lower=1000)
    df["S_SIZE"] = (np.log10(size)/np.log10(size.max())).clip(0,1)

    nc = pd.to_numeric(df["n_country"], errors="coerce").fillna(1)
    df["S_GEO"] = (np.log10(nc)/np.log10(nc.quantile(0.99))).clip(0,1)

    nh = pd.to_numeric(df["n_habitat"], errors="coerce").fillna(1)
    df["S_HAB"] = ((nh-1)/nh.quantile(0.99)).clip(0,1)

    gr = pd.to_numeric(df["annual_growth_rate"], errors="coerce").fillna(0).clip(lower=0)
    df["S_GROW"] = (gr/gr.quantile(0.99)).clip(0,1).fillna(0)

    assert abs(sum(W.values())-1.0)<1e-6, f"weights sum={sum(W.values())}"
    df["risk_score"] = (W["w_arg"]*df["S_ARG"] + W["w_mob"]*df["S_MOB"] +
                        W["w_host"]*df["S_HOST"] + W["w_rep"]*df["S_REP"] +
                        W["w_size"]*df["S_SIZE"] + W["w_geo"]*df["S_GEO"] +
                        W["w_habitat"]*df["S_HAB"] + W["w_growth"]*df["S_GROW"]).clip(0,1)

    q = G["quantiles"]
    thr = df["risk_score"].quantile(q).tolist()
    def grade(x):
        if x>=thr[2]: return "Q1"
        if x>=thr[1]: return "Q2"
        if x>=thr[0]: return "Q3"
        return "Q4"
    df["Q_grade"] = df["risk_score"].map(grade)

    out = df.sort_values("risk_score", ascending=False)
    out.to_csv(os.path.join(cfg["_out"],"psc_risk_scores.tsv"), sep="\t", index=False)
    log(f"wrote psc_risk_scores.tsv ({len(out):,} rows)")
    log("DONE.")

if __name__ == "__main__":
    main()
