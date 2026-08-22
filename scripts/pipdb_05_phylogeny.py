#!/usr/bin/env python3
"""pipdb_05_phylogeny.py — Build per-replicon core-gene alignments and ML trees."""
import argparse, os, sys, glob, subprocess
import pandas as pd

def load_config(p):
    import yaml
    with open(p) as fh: cfg = yaml.safe_load(fh)
    root = os.path.dirname(os.path.abspath(p))
    cfg["_root"]=root; cfg["_out"]=os.path.join(root,cfg["paths"]["out_dir"])
    cfg["_evo"]=os.path.join(cfg["_out"],"evolution")
    return cfg

def run(cmd, cwd=None):
    print(">>", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)

def count_fasta(fa_path):
    n = 0
    with open(fa_path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"): n += 1
    return n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--replicons", default="")
    ap.add_argument("--method", default="auto", choices=["auto","parsnp","panaroo"])
    args = ap.parse_args()
    cfg = load_config(args.config)
    P = cfg["phylogeny"]
    man = pd.read_csv(os.path.join(cfg["_evo"],"evolution_manifest.tsv"), sep="\t")
    if args.replicons:
        man = man[man["replicon"].isin(args.replicons.split(","))]
    for _, row in man.iterrows():
        rep = row["replicon"]; rd = row["dir"]
        safe = rep.replace("/","_").replace(" ","_")
        fasta_dir = os.path.join(rd,"fasta")
        multi_fa = os.path.join(fasta_dir, f"{safe}.fasta")
        gff_dir = os.path.join(rd,"gff")
        tree_dir = os.path.join(rd,"tree"); os.makedirs(tree_dir, exist_ok=True)
        gffs = sorted(glob.glob(os.path.join(gff_dir,"*","*.gff")))
        n_seqs = 0
        if os.path.exists(multi_fa): n_seqs = count_fasta(multi_fa)
        elif os.path.isdir(fasta_dir): n_seqs = len(glob.glob(os.path.join(fasta_dir,"*.fasta")))
        if n_seqs < 4:
            print(f"[skip] {rep}: no sequences found"); continue
        print(f"\n=== {rep}: {n_seqs} sequences ===")
        if args.method in ("auto","panaroo") and len(gffs) >= 4:
            pan_out = os.path.join(tree_dir,"panaroo")
            if not os.path.exists(os.path.join(pan_out,"core_gene_alignment.aln")):
                os.makedirs(pan_out, exist_ok=True)
                run(["panaroo","-i",*gffs,"-o",pan_out,"--clean-mode","strict",
                     "--core_threshold","0.95","-a","core","--aligner","mafft",
                     "--core-nucleotide","-t",str(P["threads"]),"--no_clean_edges"])
            core_gene = os.path.join(pan_out,"core_gene_alignment.aln")
            if os.path.exists(core_gene):
                prefix = os.path.join(tree_dir,f"{safe}_core")
                if not os.path.exists(prefix+".treefile"):
                    run(["iqtree2","-s",core_gene,"-m",P["iqtree_model"],"-B","1000","-T",str(P["threads"]),"--prefix",prefix])
    print("\nALL PHYLOGENIES DONE.")

if __name__ == "__main__":
    main()
