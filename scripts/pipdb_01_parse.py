#!/usr/bin/env python3
"""
pipdb_01_parse.py
=================
Integrate PIPdb PSC metadata tables into one clean master table.

Output: results/psc_master.tsv  (one row per PSC)

Usage:
  python pipdb_01_parse.py --config ../config.yaml
"""
import argparse, os, sys, json, re
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def load_config(p):
    import yaml
    with open(p) as fh:
        cfg = yaml.safe_load(fh)
    root = os.path.dirname(os.path.abspath(p))
    cfg["_root"] = root
    cfg["_db"] = os.path.join(root, cfg["paths"]["pipdb_dir"])
    cfg["_out"] = os.path.join(root, cfg["paths"]["out_dir"])
    cfg["_ref"] = os.path.join(root, cfg["paths"]["ref_dir"])
    os.makedirs(cfg["_out"], exist_ok=True)
    return cfg

def log(m): print(f"[{pd.Timestamp.now():%H:%M:%S}] {m}", flush=True)

USECOLS = [
    "id","plasmid_acc","pmlst","replicon_type","plasmid_name",
    "plasmid_seq_num","length_avg","length_min","length_max",
    "isolate_mark","isolate_mark_num","host_rank2","host_rank2_num",
    "country","country_num","collection_year","collection_year_min","collection_year_max",
    "species_name","gram_stain","genus_name","phylum_name",
    "aro_name","aro_name_num","drugclass","drugclass_num","arg_WHO","arg_WHO_num",
    "vf_name","vf_name_num","vf_category","vf_category_num",
    "gene_t4cp","gene_t4cp_num","gene_relaxase","gene_relaxase_num",
    "gene_oriT_acc","gene_oriT_num","gene_auxiliary","gene_auxiliary_num",
    "gene_IS","gene_IS_num","gene_ISfamily",
    "gene_bacmet","gene_bacmet_num",
    "annual_growth_rate","combined_risk_index",
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="../config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)
    yf = cfg["year_filter"]

    main_p = os.path.join(cfg["_db"], cfg["files"]["psc_main"])
    log("reading main PSC table ...")
    df = pd.read_csv(main_p, sep="\t", usecols=USECOLS, low_memory=False, dtype=str)
    log(f"  {len(df):,} PSCs")

    log("reading integron annotation from ARG subset ...")
    arg_p = os.path.join(cfg["_db"], cfg["files"]["psc_arg"])
    integron = pd.read_csv(arg_p, sep="\t", usecols=["plasmid_acc","gene_Integron"],
                           low_memory=False, dtype=str)
    integron["has_integron"] = integron["gene_Integron"].fillna("\\N").ne("\\N")
    df = df.merge(integron[["plasmid_acc","has_integron"]], on="plasmid_acc", how="left")
    df["has_integron"] = df["has_integron"].fillna(False)

    num_cols = ["plasmid_seq_num","length_avg","length_min","length_max",
                "isolate_mark_num","host_rank2_num","country_num",
                "collection_year_min","collection_year_max",
                "aro_name_num","drugclass_num","arg_WHO_num","vf_name_num",
                "gene_t4cp_num","gene_relaxase_num","gene_oriT_num",
                "gene_auxiliary_num","gene_IS_num","gene_bacmet_num",
                "annual_growth_rate","combined_risk_index"]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c].replace("\\N", np.nan), errors="coerce")

    df["year_min"] = df["collection_year_min"].where(
        (df["collection_year_min"]>=yf["min"]) & (df["collection_year_min"]<=yf["max"]))
    df["year_max"] = df["collection_year_max"].where(
        (df["collection_year_max"]>=yf["min"]) & (df["collection_year_max"]<=yf["max"]))
    df["year_mid"] = ((df["year_min"]+df["year_max"])/2).round(1)
    df["single_year"] = df["year_min"]==df["year_max"]

    df["n_arg"] = df["aro_name_num"].fillna(0).astype(int)
    df["n_who_arg"] = df["arg_WHO_num"].fillna(0).astype(int)
    df["n_drugclass"] = df["drugclass_num"].fillna(0).astype(int)
    df["n_vf"] = df["vf_name_num"].fillna(0).astype(int)
    df["n_is"] = df["gene_IS_num"].fillna(0).astype(int)
    df["n_metal"] = df["gene_bacmet_num"].fillna(0).astype(int)
    df["n_country"] = df["country_num"].fillna(1).astype(int)
    df["n_habitat"] = df["isolate_mark_num"].fillna(1).astype(int)

    t4cp = df["gene_t4cp_num"].fillna(0) > 0
    relax = df["gene_relaxase_num"].fillna(0) > 0
    oriT = df["gene_oriT_num"].fillna(0) > 0
    aux = df["gene_auxiliary_num"].fillna(0) > 0
    df["mobility_class"] = np.select(
        [t4cp & relax & oriT,
         t4cp & relax,
         (relax | oriT) & ~t4cp],
        ["conjugative_complete","conjugative_likely","mobilizable"],
        default="non-mobilizable")
    df["has_t4cp"] = t4cp; df["has_relaxase"] = relax
    df["has_oriT"] = oriT; df["has_auxiliary"] = aux

    def is_fam_count(x):
        if not isinstance(x,str) or x in ("\\N",""): return 0
        try: return len(json.loads(x))
        except Exception: return 0
    df["n_is_family"] = df["gene_ISfamily"].map(is_fam_count)
    df["is_density_per_kb"] = np.where(df["length_avg"]>0,
                                       df["n_is"]/(df["length_avg"]/1000), 0.0)

    df["replicon_type"] = df["replicon_type"].replace(["\\N","-",""], pd.NA)
    df["pmlst"] = df["pmlst"].replace(["\\N","-",""], pd.NA)
    rep_str = df["replicon_type"].where(df["replicon_type"].notna(), "").astype(str)
    df["replicon_primary"] = rep_str.str.split(",").str[0].replace({"":pd.NA, "nan":pd.NA, "<NA>":pd.NA})

    hr = df["host_rank2"].fillna("")
    df["host_human"] = hr.str.contains("Human", na=False)
    df["host_animal"] = hr.str.contains("Mammal|Birds|Animal", na=False, regex=True)
    im = df["isolate_mark"].fillna("")
    df["habitat_human"] = im.str.contains("Human|Feces|Blood|Urine|Sputum|Wound", na=False, regex=True)
    df["habitat_animal"] = im.str.contains("Animal", na=False)
    df["habitat_env"] = im.str.contains("Environment|Plant|Food", na=False, regex=True)

    ESCAPE = ["Enterococcus faecium","Staphylococcus aureus","Klebsiella pneumoniae",
              "Acinetobacter baumannii","Pseudomonas aeruginosa","Enterobacter"]
    WHO_CRITICAL = ["Escherichia coli","Salmonella enterica","Shigella",
                    "Citrobacter","Serratia","Proteus","Morganella","Providencia",
                    "Klebsiella aerogenes","Klebsiella oxytoca","Enterobacter hormaechei",
                    "Enterobacter cloacae","Campylobacter"]
    def host_class(s):
        if not isinstance(s,str): return "unknown"
        for e in ESCAPE:
            if e in s: return "escape"
        for w in WHO_CRITICAL:
            if w in s: return "who_critical"
        return "pathogen_other"
    df["host_class"] = df["species_name"].map(host_class)

    out_cols = ["id","plasmid_acc","plasmid_name","replicon_type","replicon_primary",
                "pmlst","length_avg","length_min","length_max","plasmid_seq_num",
                "species_name","genus_name","phylum_name","gram_stain","host_class",
                "host_rank2","host_human","host_animal",
                "isolate_mark","n_habitat","habitat_human","habitat_animal","habitat_env",
                "country","n_country",
                "year_min","year_max","year_mid","single_year",
                "aro_name","n_arg","n_who_arg","drugclass","n_drugclass",
                "vf_name","vf_category","n_vf",
                "gene_t4cp","gene_relaxase","gene_oriT_acc","gene_auxiliary",
                "has_t4cp","has_relaxase","has_oriT","has_auxiliary","has_integron",
                "mobility_class",
                "gene_IS","n_is","n_is_family","is_density_per_kb",
                "gene_bacmet","n_metal",
                "annual_growth_rate","combined_risk_index"]
    df[out_cols].to_csv(os.path.join(cfg["_out"],"psc_master.tsv"), sep="\t", index=False)
    log(f"wrote psc_master.tsv: {len(df):,} rows, {len(out_cols)} columns")
    log("DONE.")

if __name__ == "__main__":
    main()
