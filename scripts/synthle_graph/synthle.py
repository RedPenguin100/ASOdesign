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
MIN_DEG = 1          # discard tiny nodes
MAX_DEG_PERCENTILE = 99.0  # discard huge hubs
LOW_THROUGHPUT = True

# ==============


def load_biogrid_tab2(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # human-human only
    df = df[(df["Organism Interactor A"] == 9606) & (df["Organism Interactor B"] == 9606)]
    df["Experimental System Type"] = df["Experimental System Type"].str.lower()
    df["Throughput"] = df["Throughput"].fillna("").str.title()
    return df


def filter_df(df, type="genetic"):
    
    if type == "physical":
        df = df[df["Experimental System Type"] == "physical"]
    elif type == "genetic":
        df = df[df["Experimental System Type"] == "genetic"]
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
    TYPE = "genetic"
    OUT_CSV = f"tp53_{TYPE}_results.csv"

    df0 = load_biogrid_tab2(BIOGRID_FILE)
    df1 = filter_df(df0, TYPE)
    edf = collapse_edges(df1)
    partners = build_partner_sets(edf)
    
    print(df1.head())
    # Count how many candidates have an edge with TP53 in "genetic" type
    TP53 = "TP53"
    if TP53 in partners:
        tp53_partners = partners[TP53]
        print(f"Number of candidates with an edge to TP53 in 'genetic' type: {len(tp53_partners)}")
    else:
        print("TP53 not found in partners for 'genetic' type.")

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
    # Add all candidate proteins in the first column and a column indicating if they are TP53 neighbors
    all_candidates = list(universe)
    # For each candidate, check if it is a TP53 neighbor
    is_tp53_neighbor = [1 if c in tp53_partners else 0 for c in all_candidates]
    # Create a DataFrame for all candidates and their TP53 neighbor status
    candidates_df = pd.DataFrame({
        "CandidateProtein": all_candidates,
        "IsTP53Neighbor": is_tp53_neighbor
    })
    # Merge the results DataFrame with the candidates DataFrame
    res = candidates_df.merge(res, left_on="CandidateProtein", right_on="Protein", how="left")
    res = res.drop(columns=["Protein"])

    res = res.sort_values(["HypergeomP", "Jaccard", "Degree"],
                          ascending=[True, False, False]).reset_index(drop=True)

    print(f"TP53 degree={target_deg}, candidates={len(res)}")
    filtered_res = res
    # genes = [
    # "BGN", "TP53BP1", "BCAM", "MBOAT1", "USP28", "GABRR3", "DDX3Y", "CD38",
    # "QPCT", "CDKN1A", "CDC14A", "UBQLN2", "FEM1B", "GPR31", "RAB5B", "JAZF1",
    # "PNPLA6", "RGS12", "RARA", "LANCL2", "FADD", "APTX", "ZBTB4", "ZRSR2",
    # "GTPBP6", "RAB6A", "WTH3DI", "NFE2L2", "AHCTF1", "MAPKAPK2", "CTNNB1",
    # "MET", "SGK2", "RB1", "CSNK1E", "PAK3"
    # ]
    # Display only the genes from the list
    # filtered_res = res[res["CandidateProtein"].isin(genes)].reset_index(drop=True)
    # quick view
    filtered_res = filtered_res.sort_values("HypergeomP", ascending=True).reset_index(drop=True)
    print(filtered_res.head(TOPK).to_string(index=False))
    
    filtered_res.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {TYPE} results to {OUT_CSV}")
    print(filtered_res.head(TOPK).to_string(index=False))


if __name__ == "__main__":
    main()
