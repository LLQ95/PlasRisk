#!/usr/bin/env python3
"""pipdb_03_evolution.py — Prepare tip-dated datasets for per-replicon phylogeny."""
import argparse, os, sys, subprocess, time
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def load_config(p):
    import yaml
    with open(p) as fh: cfg = yaml.safe_load(fh)
    root = os.path.dirname(os.path.abspath(p))
    cfg["_root"]=root; cfg["_db"]=os.path.join(root,cfg["paths"]["pipdb_dir"])
    cfg["_out"]=os.path.join(root,cfg["paths"]["out_dir"])
    cfg["_fasta"]=os.path.join(root,cfg["paths"]["fasta_dir"])
    os.makedirs(cfg["_out"], exist_ok=True); os.makedirs(cfg["_fasta"], exist_ok=True)
    return cfg

def log(m): print(f"[{pd.Timestamp.now():%H:%M:%S}] {m}", flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="../config.yaml")
    ap.add_argument("--replicons", default="")
    ap.add_argument("--download", action="store_true")
    args = ap.parse_args()
    cfg = load_config(args.config)
    P = cfg["phylogeny"]; yf = cfg["year_filter"]
    risk = pd.read_csv(os.path.join(cfg["_out"],"psc_risk_scores.tsv"), sep="\t",
                       usecols=["plasmid_acc","Q_grade","replicon_primary","risk_score"], low_memory=False)
    risk_map = dict(zip(risk["plasmid_acc"], risk["Q_grade"]))
    cpath = os.path.join(cfg["_db"], cfg["files"]["contig"])
    c = pd.read_csv(cpath, sep="\t", low_memory=False, dtype=str,
                    usecols=["plasmid_seq_id","plasmid_acc","replicon_type","pmlst",
                             "seq_id","assembly_acc","length","country","collection_year",
                             "isolate_mark","host_rank2","species_name","genus_name"])
    c["length"] = pd.to_numeric(c["length"], errors="coerce")
    c["year"] = pd.to_numeric(c["collection_year"], errors="coerce")
    c = c[(c["year"]>=yf["min"])&(c["year"]<=yf["max"])&(c["length"]>=P["min_length_bp"])]
    c["Q_grade"] = c["plasmid_acc"].map(risk_map).fillna("Q4")
    c["replicon_primary"] = c["replicon_type"].replace(["\\N","-",""],pd.NA).astype(str).str.split(",").str[0]
    c.loc[c["replicon_primary"].isin(["nan","None","<NA>",""]), "replicon_primary"] = pd.NA
    c = c[c["replicon_primary"].notna()]
    reps = [r.strip() for r in args.replicons.split(",")] if args.replicons else \
        (c["replicon_primary"].value_counts().loc[lambda s: s>=P["min_plasmids_per_replicon"]].index.tolist())
    evo_dir = os.path.join(cfg["_out"],"evolution"); os.makedirs(evo_dir, exist_ok=True)
    manifest = []
    for rep in reps:
        sub = c[c["replicon_primary"]==rep].copy()
        if len(sub) < P["min_plasmids_per_replicon"]: continue
        if len(sub) > P["max_per_replicon"]:
            sub["decade"] = (sub["year"]//10*10).astype(int)
            frac = P["max_per_replicon"]/len(sub)
            sub = (sub.groupby(["decade","Q_grade"], group_keys=False)
                      .sample(frac=min(1.0, frac*1.5), random_state=42))
            sub = sub.reset_index(drop=True).drop_duplicates("plasmid_seq_id").head(P["max_per_replicon"])
        rd = os.path.join(evo_dir, rep.replace("/","_").replace(" ","_"))
        os.makedirs(rd, exist_ok=True)
        out = sub[["plasmid_seq_id","seq_id","assembly_acc","plasmid_acc",
                   "species_name","genus_name","country","isolate_mark","host_rank2",
                   "year","length","Q_grade"]].copy()
        out.to_csv(os.path.join(rd,"tip_dates.tsv"), sep="\t", index=False)
        manifest.append({"replicon":rep,"n":len(out),"year_min":int(out["year"].min()),
                         "year_max":int(out["year"].max()),"dir":rd})
    pd.DataFrame(manifest).to_csv(os.path.join(evo_dir,"evolution_manifest.tsv"), sep="\t", index=False)
    log("DONE.")

if __name__ == "__main__":
    main()
