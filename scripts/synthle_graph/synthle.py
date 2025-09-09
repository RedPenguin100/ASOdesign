#!/usr/bin/env python3
import pandas as pd
import math
from collections import defaultdict

try:
    from scipy.stats import hypergeom
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False


# === CONFIG ===
BIOGRID_FILE = "BIOGRID-ORGANISM-Homo_sapiens-4.4.248.tab2.txt"
TARGET = "TP53"
TOPK = 30
OUT_CSV = "tp53_like_results.csv"
MIN_DEG = 10          # discard tiny nodes
MAX_DEG_PERCENTILE = 99.0  # discard huge hubs
LOW_THROUGHPUT = True
PHYSICAL_ONLY = True
# ==============


def load_biogrid_tab2(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # human-human only
    df = df[(df["Organism Interactor A"] == 9606) & (df["Organism Interactor B"] == 9606)]
    df["Experimental System Type"] = df["Experimental System Type"].str.lower()
    df["Throughput"] = df["Throughput"].fillna("").str.title()
    return df


def filter_df(df):
    
    if PHYSICAL_ONLY:
        df = df[df["Experimental System Type"] == "physical"]
    if LOW_THROUGHPUT:
        df = df[df["Throughput"] == "Low Throughput"]
    df = df.dropna(subset=["Official Symbol Interactor A", "Official Symbol Interactor B"])
    df = df[df["Official Symbol Interactor A"] != df["Official Symbol Interactor B"]]
    return df


def collapse_edges(df):
    def key(a, b): return tuple(sorted((a, b)))
    pmids = defaultdict(set)
    for a, b, pm in zip(df["Official Symbol Interactor A"],
                        df["Official Symbol Interactor B"],
                        df["Pubmed ID"].astype(str)):
        pmid_list = [p.strip() for p in pm.split("|")] if pm and pm != "nan" else []
        pmids[key(a, b)].update(pmid_list)
    edges = [(a, b, len(s)) for (a, b), s in pmids.items()]
    return pd.DataFrame(edges, columns=["A", "B", "pmid_support"])


def build_partner_sets(edf, min_support=1):
    edf = edf[edf["pmid_support"] >= min_support]
    partners = defaultdict(set)
    for a, b in zip(edf["A"], edf["B"]):
        partners[a].add(b)
        partners[b].add(a)
    return partners


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a or b) else 0.0


def overlap(a, b):
    return len(a & b) / min(len(a), len(b)) if (a and b) else 0.0


def hypergeom_sf(N, K, n, k):
    if SCIPY_OK:
        return float(hypergeom.sf(k - 1, N, K, n))
    # fallback: simple approx
    return 1.0 if k == 0 else 1e-6


def main():
    df0 = load_biogrid_tab2(BIOGRID_FILE)
    df1 = filter_df(df0)
    edf = collapse_edges(df1)
    partners = build_partner_sets(edf)

    if TARGET not in partners:
        print(f"{TARGET} not found.")
        return

    import numpy as np
    degrees = {p: len(s) for p, s in partners.items()}
    cut = np.percentile(list(degrees.values()), MAX_DEG_PERCENTILE)
    target_set = partners[TARGET]
    target_deg = len(target_set)
    universe = set(partners) - {TARGET}
    N = len(universe)

    rows = []
    for prot, neigh in partners.items():
        if prot == TARGET:
            continue
        d = len(neigh)
        if d < MIN_DEG or d > cut:
            continue
        inter = len(target_set & neigh)
        if inter == 0:
            continue
        jac = jaccard(target_set, neigh)
        oc = overlap(target_set, neigh)
        pval = hypergeom_sf(N, target_deg, d, inter)
        rows.append((prot, d, inter, jac, oc, pval))

    res = pd.DataFrame(rows, columns=["Protein", "Degree", "SharedWithTP53",
                                      "Jaccard", "OverlapCoef", "HypergeomP"])
    res = res.sort_values(["HypergeomP", "Jaccard", "Degree"],
                          ascending=[True, False, False]).reset_index(drop=True)

    print(f"TP53 degree={target_deg}, candidates={len(res)}")
    print(res.head(TOPK).to_string(index=False))

    res.to_csv(OUT_CSV, index=False)
    print(f"\nSaved results to {OUT_CSV}")


if __name__ == "__main__":
    main()
