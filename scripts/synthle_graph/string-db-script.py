import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import hypergeom

def tp53_like_single_channel(
    edge_type: str,
    info_df: pd.DataFrame,
    df: pd.DataFrame,
    target: str = "TP53",
    min_combined: int = 0,
    min_deg: int = 1,
    max_deg_percentile: float = 99.0,
    topk: int = 30,
) -> pd.DataFrame:
    
    et = edge_type
    allowed = {
        "neighborhood","fusion","cooccurence","coexpression","experimental","database","textmining, combined_score"
    }
    if et not in allowed:
        raise ValueError(f"edge_type must be one of {sorted(allowed)}; got '{edge_type}'")

    # --- score + single-channel filter ---
    work = df.copy()
    work["combined_score"] = pd.to_numeric(work["combined_score"], errors="coerce")
    work = work.dropna(subset=["combined_score"])
    
    if et not in work.columns:
        raise ValueError(f"Column '{et}' not found in df; available: {sorted(work.columns)}")
    
    work[et] = pd.to_numeric(work[et], errors="coerce").fillna(0)

    work = work[(work["combined_score"] >= min_combined) & (work[et] > 0)]
    
    if work.empty:
        raise ValueError(f"No edges remain for edge_type='{et}' at combined_score >= {min_combined}")

    # --- map STRING IDs -> gene symbols ---
    id_to_gene = dict(zip(info_df["protein"], info_df["gene"]))
    g = work.assign(A=work["protein1"].map(id_to_gene),
                    B=work["protein2"].map(id_to_gene)).dropna(subset=["A","B"])
    g = g[g["A"] != g["B"]]

    # --- undirected canonicalization + collapse duplicates (max over numeric cols) ---
    swap = g["A"] > g["B"]
    g.loc[swap, ["A","B"]] = g.loc[swap, ["B","A"]].values

    agg = {"combined_score": "max", et: "max"}
    for c in g.columns:
        if c in ("A","B","protein1","protein2", et, "combined_score"): 
            continue
        if pd.api.types.is_numeric_dtype(g[c]):
            agg[c] = "max"
    g = g.groupby(["A","B"], as_index=False).agg(agg)

    # --- neighbor sets ---
    neighbors = defaultdict(set)
    for a,b in zip(g["A"], g["B"]):
        neighbors[a].add(b); neighbors[b].add(a)
    if target not in neighbors:
        raise ValueError(f"{target} not present after filtering for edge_type='{et}'.")

    # --- similarity metrics ---
    def jacc(a,b): return len(a & b) / len(a | b) if (a or b) else 0.0
    def overlap(a,b): return len(a & b) / min(len(a), len(b)) if (a and b) else 0.0

  

    def hypergeom_sf(N, K, n, k):
        # Returns the survival function (P[X >= k]) for the hypergeometric distribution
        # N: population size, K: number of "successes" in population,
        # n: number of draws, k: observed successes
        return float(hypergeom.sf(k-1, N, K, n))

    # --- rank TP53-like genes ---
    deg = {x: len(s) for x,s in neighbors.items()}
    cut = np.percentile(list(deg.values()), max_deg_percentile)

    T = neighbors[target]; K = len(T)
    U = set(neighbors) - {target}; N = len(U)

    rows = []
    for gene, nbrs in neighbors.items():
        if gene == target: continue
        d = len(nbrs)
        if d < min_deg or d > cut: continue
        inter = len(T & nbrs)
        if inter == 0: continue
        rows.append([gene, d, inter, jacc(T, nbrs), overlap(T, nbrs), hypergeom_sf(N, K, d, inter)])

    res = (pd.DataFrame(rows, columns=["Gene","Degree","SharedWithTP53","Jaccard","OverlapCoef","HypergeomP"])
            .sort_values(["HypergeomP","Jaccard","Degree"], ascending=[True,False,False])
            .reset_index(drop=True))
    
    # Add all candidate genes in the first column and a column indicating if they are TP53 neighbors
    all_candidates = list(U)
    is_tp53_neighbor = [1 if c in T else 0 for c in all_candidates]
    candidates_df = pd.DataFrame({
        "CandidateGene": all_candidates,
        "IsTP53Neighbor": is_tp53_neighbor
    })
    # Merge the results DataFrame with the candidates DataFrame
    res = candidates_df.merge(res, left_on="CandidateGene", right_on="Gene", how="left")
    res = res.drop(columns=["Gene"])
    
    filtered_res = res
    # genes = [
    # "BGN", "TP53BP1", "BCAM", "MBOAT1", "USP28", "GABRR3", "DDX3Y", "CD38",
    # "QPCT", "CDKN1A", "CDC14A", "UBQLN2", "FEM1B", "GPR31", "RAB5B", "JAZF1",
    # "PNPLA6", "RGS12", "RARA", "LANCL2", "FADD", "APTX", "ZBTB4", "ZRSR2",
    # "GTPBP6", "RAB6A", "WTH3DI", "NFE2L2", "AHCTF1", "MAPKAPK2", "CTNNB1",
    # "MET", "SGK2", "RB1", "CSNK1E", "PAK3"
    # ]
    # # Display only the genes from the list
    # filtered_res = res[res["CandidateGene"].isin(genes)].reset_index(drop=True)
    
    # quick view
    filtered_res = filtered_res.sort_values("HypergeomP", ascending=True).reset_index(drop=True)
    print(f"[edge_type={edge_type}] TP53 degree={K}, candidates={len(filtered_res)}, edges_kept={len(g)} (min_combined={min_combined})")
    print(filtered_res.head(topk).to_string(index=False))
    
    filtered_res.to_csv(OUT_CSV, index=False)
    print(f"\nSaved {edge_type} results to {OUT_CSV}")
    
    return res

# --- example usage (assumes you already ran the two lines you posted) ---
info_df = pd.read_csv("string_info.csv")
df = pd.read_csv("string_detailed.csv")

allowed = {"neighborhood","fusion","cooccurence","coexpression","experimental","database","textmining, combined_score"}

for edge_type in allowed:
    print(f"\n--- Running tp53_like_single_channel for edge_type: {edge_type} ---")
    try:
        OUT_CSV = f"tp53_{edge_type}_results.csv"
        res = tp53_like_single_channel(edge_type, info_df, df)
        
    except Exception as e:
        print(f"Error for edge_type '{edge_type}': {e}")
        continue


