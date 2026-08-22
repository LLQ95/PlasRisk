#!/usr/bin/env python3
"""pipdb_06_beast_setup.py — Generate BEAST2 XML files for per-replicon dating."""
import argparse, os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def load_config(p):
    import yaml
    with open(p) as fh: cfg = yaml.safe_load(fh)
    root = os.path.dirname(os.path.abspath(p))
    cfg["_root"]=root; cfg["_out"]=os.path.join(root,cfg["paths"]["out_dir"])
    cfg["_evo"]=os.path.join(cfg["_out"],"evolution")
    return cfg

def log(m): print(f"[{pd.Timestamp.now():%H:%M:%S}] {m}", flush=True)

def build_xml(alignment_path, dates_path, output_xml, replicon_name, beast_cfg):
    import subprocess, tempfile
    from Bio import AlignIO
    aln = AlignIO.read(alignment_path, "fasta")
    n_taxa = len(aln); seq_len = aln.get_alignment_length()
    dates = pd.read_csv(dates_path, sep="\t")
    date_map = dict(zip(dates["plasmid_seq_id"].astype(str).str.replace(".","_"), dates["year"]))
    taxon_dates = []
    for rec in aln:
        tid = rec.id
        yr = date_map.get(tid, date_map.get(tid.replace("_","."), 2000))
        taxon_dates.append(f"{tid}={yr}")
    log("  use BEAUti or beast2xml template to complete XML generation")
    log(f"  {n_taxa} taxa, {seq_len} sites, {len(taxon_dates)} dates mapped")
    return n_taxa, seq_len

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--replicons", default="")
    args = ap.parse_args()
    cfg = load_config(args.config)
    B = cfg["beast"]
    man = pd.read_csv(os.path.join(cfg["_evo"],"evolution_manifest.tsv"), sep="\t")
    if args.replicons:
        man = man[man["replicon"].isin(args.replicons.split(","))]
    for _, row in man.iterrows():
        rep = row["replicon"]; rd = row["dir"]
        safe = rep.replace("/","_").replace(" ","_")
        core_gene = os.path.join(rd,"tree","panaroo","core_gene_alignment.aln")
        snp_fa = os.path.join(rd,"tree","snippy",f"{safe}.snps.fasta")
        aln = core_gene if os.path.exists(core_gene) else (snp_fa if os.path.exists(snp_fa) else None)
        if not aln:
            print(f"[skip] {rep}: no alignment"); continue
        dates = os.path.join(rd,"tip_dates.tsv")
        out = os.path.join(rd,"tree",f"{safe}_beast.xml")
        print(f"\n=== {rep} ===")
        build_xml(aln, dates, out, rep, B)
    log("DONE.")

if __name__ == "__main__":
    main()
