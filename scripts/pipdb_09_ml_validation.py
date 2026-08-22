#!/usr/bin/env python3
"""pipdb_09_ml_validation.py — Random Forest + SHAP validation of risk components."""
import argparse, os
import pandas as pd
import numpy as np

def load_config(p):
    import yaml
    with open(p) as fh: cfg = yaml.safe_load(fh)
    root = os.path.dirname(os.path.abspath(p))
    cfg["_out"]=os.path.join(root,cfg["paths"]["out_dir"])
    return cfg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="../config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import train_test_split
    df = pd.read_csv(os.path.join(cfg["_out"],"psc_risk_scores.tsv"), sep="\t", low_memory=False)
    feats = ["S_ARG","S_MOB","S_HOST","S_REP","S_SIZE","S_GEO","S_HAB","S_GROW",
             "n_arg","n_who_arg","n_high_risk_arg","n_vf","n_is","n_metal",
             "n_country","n_habitat","length_avg","annual_growth_rate","has_integron"]
    for f in feats:
        if f not in df.columns: df[f]=0
        df[f] = pd.to_numeric(df[f], errors="coerce").fillna(0)
    X = df[feats]
    Xtr,Xte,ytr,yte = train_test_split(X, df["risk_score"], test_size=0.2, random_state=42)
    rf = RandomForestRegressor(n_estimators=200, n_jobs=-1, random_state=42, max_depth=12)
    rf.fit(Xtr,ytr)
    print(f"R^2: {rf.score(Xte,yte):.3f}")
    imp = pd.DataFrame({"feature":feats,"importance":rf.feature_importances_}).sort_values("importance",ascending=False)
    imp.to_csv(os.path.join(cfg["_out"],"ml_feature_importance.tsv"), sep="\t", index=False)
    print(imp.to_string(index=False))

if __name__ == "__main__":
    main()
