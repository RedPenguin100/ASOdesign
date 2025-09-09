#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rank synthetic-lethal partners for a target gene (default: TP53) from a LOCAL TSV only.

Input file (default):
  /home/espl/ASOdesign/scripts/gene_sl_gene.tsv

Expected columns (case-insensitive; flexible):
  x:START_ID, x_type, x_name, x_source, y:END_ID, y_type, y_name, y_source,
  relation, :TYPE, rel_source, edge_index, cell_line, pubmed_id, cancer

Key upgrades vs previous version
--------------------------------
1) Weighted evidence: rel_source/assay mapped to confidence weights (CRISPR>RNAi>Computational).
2) Score incorporates weighted evidence and unique PMIDs/cell lines.
3) CLI:
   --target/-t, --input/-i, --out-prefix/-o
   --top-k/-k (preview lines), --min-evidence (filter weak partners),
   --explain PARTNER (print PMIDs, sources, cell-line breakdown).
4) Stable tie-breaking.

Usage
-----
  python synthle.py
  python synthle.py -t TP53 -i /home/espl/ASOdesign/scripts/gene_sl_gene.tsv -o tp53 --top-k 20
  python synthle.py --explain MTOR

Outputs
-------
  <prefix>_partners_ranked.csv
  <prefix>_partners_topK.csv
  <prefix>_partners_by_cell_line.csv   (if cell_line present)
"""

import argparse
import math
from pathlib import Path
from typing import Dict, Tuple, List

import pandas as pd


# ---------- Defaults ----------
DEFAULT_INPUT = "/home/espl/ASOdesign/scripts/gene_sl_gene.tsv"
DEFAULT_TARGET = "TP53"
DEFAULT_PREFIX = "tp53"
DEFAULT_TOPK = 30
DEFAULT_MIN_EVIDENCE = 1  # drop partners with fewer total rows than this


# ---------- Evidence weight map ----------
# Map rel_source / assay strings to weights. You can tune freely.
WEIGHT_MAP: Dict[str, float] = {
    "crispr": 1.0,
    "crispri": 1.0,
    "crisper": 1.0,  # common typo guard
    "rnai": 0.7,
    "sirna": 0.7,
    "shrna": 0.7,
    "genomernai": 0.7,
    "osprey": 0.8,
    "y2h": 0.6,
    "computational": 0.25,
    "computational prediction": 0.25,
    "prediction": 0.25,
    # fallbacks below
}

DEFAULT_WEIGHT = 0.5  # used if no mapping matched


def weight_for_source(s: str) -> float:
    if not isinstance(s, str) or not s:
        return DEFAULT_WEIGHT
    l = s.lower()
    # try direct match
    if l in WEIGHT_MAP:
        return WEIGHT_MAP[l]
    # substring matches
    for k, w in WEIGHT_MAP.items():
        if k in l:
            return w
    return DEFAULT_WEIGHT


# ---------- Helpers ----------
def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]
    return df


def pick_column(df: pd.DataFrame, candidates) -> str:
    lower_map = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def load_edges(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, encoding="utf-8")
    df = normalize_columns(df)

    x_name = pick_column(df, ["x_name", "x name", "genea", "gene a", "a", "source_gene"])
    y_name = pick_column(df, ["y_name", "y name", "geneb", "gene b", "b", "target_gene"])
    relation = pick_column(df, ["relation"])
    type_col = pick_column(df, [":TYPE", "type", "edge_type"])
    rel_source = pick_column(df, ["rel_source", "source", "method", "assay"])
    cell_line = pick_column(df, ["cell_line", "cell line", "cellline"])
    pubmed = pick_column(df, ["pubmed_id", "pubmed", "pmid"])

    df.attrs["col_x_name"] = x_name
    df.attrs["col_y_name"] = y_name
    df.attrs["col_relation"] = relation
    df.attrs["col_type"] = type_col
    df.attrs["col_rel_source"] = rel_source
    df.attrs["col_cell_line"] = cell_line
    df.attrs["col_pubmed"] = pubmed

    if x_name is None or y_name is None:
        raise ValueError(
            "Could not find gene symbol columns (x_name / y_name). "
            f"Columns present: {list(df.columns)}"
        )

    # Normalize
    df[x_name] = df[x_name].astype(str).str.strip().str.upper()
    df[y_name] = df[y_name].astype(str).str.strip().str.upper()
    for col in [relation, type_col, rel_source, cell_line, pubmed]:
        if col and col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    return df


def is_synthetic_lethal_row(row, relation_col: str, type_col: str) -> bool:
    rel_ok = False
    typ_ok = False
    if relation_col and relation_col in row and isinstance(row[relation_col], str):
        rel_ok = "sl" in row[relation_col].lower()
    if type_col and type_col in row and isinstance(row[type_col], str):
        typ_ok = "gene_sl_gene" in row[type_col].lower()
    return rel_ok or typ_ok or (relation_col is None and type_col is None)


def filter_to_SL(df: pd.DataFrame) -> pd.DataFrame:
    relation_col = df.attrs.get("col_relation")
    type_col = df.attrs.get("col_type")
    mask = df.apply(lambda r: is_synthetic_lethal_row(r, relation_col, type_col), axis=1)
    return df.loc[mask].copy()


def expand_pmids(series: pd.Series) -> pd.Series:
    """Split multi-PMID fields into individual IDs; return flattened series."""
    if series.empty:
        return series
    s = series.dropna().astype(str).str.strip()
    if s.empty:
        return s
    return s.str.split(r"[;,| ]+").explode().str.strip().replace({"": None}).dropna()


def partners_for_target(df: pd.DataFrame, target: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    x = df.attrs["col_x_name"]
    y = df.attrs["col_y_name"]
    rel_source = df.attrs.get("col_rel_source")
    cell_line = df.attrs.get("col_cell_line")
    pubmed = df.attrs.get("col_pubmed")

    use_cols = [x, y]
    if rel_source: use_cols.append(rel_source)
    if cell_line:  use_cols.append(cell_line)
    if pubmed:     use_cols.append(pubmed)

    sub = df.loc[(df[x] == target) | (df[y] == target), use_cols].copy()
    if sub.empty:
        return (
            pd.DataFrame(columns=["partner", "weighted_evidence", "evidence_count",
                                  "n_cell_lines", "n_sources", "n_pubmed", "score"]),
            pd.DataFrame()
        )

    # Identify partner
    sub["partner"] = sub.apply(lambda r: r[y] if r[x] == target else r[x], axis=1)

    # Clean columns
    if cell_line and cell_line in sub.columns:
        sub[cell_line] = sub[cell_line].replace({'""': "", '"': ""}).fillna("").astype(str).str.strip()
    if rel_source and rel_source in sub.columns:
        sub[rel_source] = sub[rel_source].fillna("").astype(str).str.strip()
    if pubmed and pubmed in sub.columns:
        sub[pubmed] = sub[pubmed].fillna("").astype(str).str.strip()

    # Per-row weight
    if rel_source and rel_source in sub.columns:
        sub["_weight"] = sub[rel_source].map(weight_for_source)
    else:
        sub["_weight"] = DEFAULT_WEIGHT

    # Aggregations
    grp = sub.groupby("partner", dropna=False)
    evidence_count = grp.size().rename("evidence_count")
    weighted_evidence = grp["_weight"].sum().rename("weighted_evidence")

    n_sources = pd.Series(0, index=grp.size().index, name="n_sources")
    if rel_source and rel_source in sub.columns:
        n_sources = grp[rel_source].apply(lambda s: s.replace("", pd.NA).dropna().str.lower().nunique()).rename("n_sources")

    n_cell_lines = pd.Series(0, index=grp.size().index, name="n_cell_lines")
    if cell_line and cell_line in sub.columns:
        n_cell_lines = grp[cell_line].apply(lambda s: s.replace("", pd.NA).dropna().nunique()).rename("n_cell_lines")

    n_pubmed = pd.Series(0, index=grp.size().index, name="n_pubmed")
    if pubmed and pubmed in sub.columns:
        # explode PMIDs then count unique per partner
        pm = sub.loc[sub[pubmed].str.len() > 0, ["partner", pubmed]].copy()
        if not pm.empty:
            pm = pm.assign(_pmid=expand_pmids(pm[pubmed]))
            pm = pm[pd.notna(pm["_pmid"])]
            if not pm.empty:
                n_pubmed = pm.groupby("partner")["_pmid"].nunique().rename("n_pubmed")

    out = pd.concat([weighted_evidence, evidence_count, n_cell_lines, n_sources, n_pubmed], axis=1).fillna(0)
    out.reset_index(inplace=True)

    # Score
    out["score"] = (
        out["weighted_evidence"].apply(lambda v: math.log1p(v))
        + 0.5 * out["n_cell_lines"].apply(lambda v: math.log1p(v))
        + 0.25 * out["n_sources"]
        + 0.5 * out["n_pubmed"].apply(lambda v: math.log1p(v))
    )

    # Stable sorting / tie-breakers
    out = out.sort_values(
        by=["score", "weighted_evidence", "evidence_count", "partner"],
        ascending=[False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    return out, sub  # return sub for explain/pivot


def pivot_by_cell_line(sub_edges: pd.DataFrame, target: str, df_attrs) -> pd.DataFrame:
    x = df_attrs["col_x_name"]
    y = df_attrs["col_y_name"]
    cell_line = df_attrs.get("col_cell_line")
    if not cell_line or cell_line not in sub_edges.columns:
        return pd.DataFrame()

    tmp = sub_edges[[x, y, cell_line]].copy()
    tmp["partner"] = tmp.apply(lambda r: r[y] if r[x] == target else r[x], axis=1)
    tmp[cell_line] = tmp[cell_line].replace({'""': "", '"': ""}).fillna("").astype(str).str.strip()
    tmp = tmp[tmp[cell_line].str.len() > 0]
    if tmp.empty:
        return pd.DataFrame()

    ct = tmp.groupby(["partner", cell_line]).size().rename("count").reset_index()
    pivot = ct.pivot(index="partner", columns=cell_line, values="count").fillna(0).astype(int)
    pivot = pivot.sort_index()
    return pivot


def explain_partner(partner: str, sub_edges: pd.DataFrame, df_attrs) -> str:
    """Return a multi-line string with details for a given partner (PMIDs, sources, cell lines)."""
    partner = partner.strip().upper()
    x = df_attrs["col_x_name"]
    y = df_attrs["col_y_name"]
    rel_source = df_attrs.get("col_rel_source")
    cell_line = df_attrs.get("col_cell_line")
    pubmed = df_attrs.get("col_pubmed")

    rows = sub_edges[(sub_edges[x] == partner) | (sub_edges[y] == partner)].copy()
    if rows.empty:
        # when partner is only on the opposite side against TARGET, recompute partner like before
        rows = sub_edges.copy()
        rows["partner"] = rows.apply(lambda r: r[y] if r[x] == partner else (r[x] if r[y] == partner else None), axis=1)
        rows = rows[rows["partner"] == partner]

    if rows.empty:
        return f"No rows found for partner {partner}."

    # Aggregate details
    sources = set()
    cells = set()
    pmids = set()

    if rel_source and rel_source in rows.columns:
        sources = set([s.strip() for s in rows[rel_source].dropna().astype(str) if s.strip()])

    if cell_line and cell_line in rows.columns:
        cells = set([s.strip() for s in rows[cell_line].dropna().astype(str) if s.strip()])

    if pubmed and pubmed in rows.columns:
        flat_pmids = expand_pmids(rows[pubmed])
        pmids = set(flat_pmids.astype(str)) if not flat_pmids.empty else set()

    msg = [
        f"Partner: {partner}",
        f"  Evidence rows: {len(rows)}",
        f"  Sources: {', '.join(sorted(sources)) if sources else '(none)'}",
        f"  Cell lines: {', '.join(sorted(cells)) if cells else '(none)'}",
        f"  PMIDs: {', '.join(sorted(pmids)) if pmids else '(none)'}",
    ]
    return "\n".join(msg)


# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Rank synthetic-lethal partners from a local TSV.")
    ap.add_argument("--input", "-i", default=DEFAULT_INPUT, help="Path to gene_sl_gene.tsv")
    ap.add_argument("--target", "-t", default=DEFAULT_TARGET, help="Target gene symbol (default: TP53)")
    ap.add_argument("--out-prefix", "-o", default=DEFAULT_PREFIX, help="Prefix for output files")
    ap.add_argument("--top-k", "-k", type=int, default=DEFAULT_TOPK, help="Rows to include in the *_topK.csv")
    ap.add_argument("--min-evidence", type=int, default=DEFAULT_MIN_EVIDENCE,
                    help="Drop partners with fewer than this many evidence rows")
    ap.add_argument("--explain", metavar="GENE", help="Print detailed evidence for a specific partner and exit")
    args = ap.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    target = args.target.strip().upper()
    prefix = args.out_prefix.strip()

    # Load & filter to SL
    edges = load_edges(str(input_path))
    edges_sl = filter_to_SL(edges)

    # Rank
    ranked, sub_edges = partners_for_target(edges_sl, target=target)

    if ranked.empty:
        print(f"No SL partners found for {target}.")
        return

    # Apply min evidence filter (on raw count, not weighted)
    if args.min_evidence > 1:
        before = len(ranked)
        ranked = ranked[ranked["evidence_count"] >= args.min_evidence].reset_index(drop=True)
        print(f"Filtered partners by min_evidence={args.min_evidence}: {before} -> {len(ranked)}")

    # Save
    out_full = Path(f"{prefix}_partners_ranked.csv")
    out_top = Path(f"{prefix}_partners_top{args.top_k}.csv")
    ranked.to_csv(out_full, index=False)
    ranked.head(args.top_k).to_csv(out_top, index=False)

    # Pivot by cell line
    pivot = pivot_by_cell_line(sub_edges, target=target, df_attrs=edges.attrs)
    out_by_cl = Path(f"{prefix}_partners_by_cell_line.csv")
    if not pivot.empty:
        pivot.to_csv(out_by_cl)

    # Preview
    print("Saved:")
    print(" ", out_full.resolve())
    print(" ", out_top.resolve())
    if not pivot.empty:
        print(" ", out_by_cl.resolve(), "(partner × cell_line counts)")
    else:
        print("  (No cell_line info to pivot)")

    # Optional explain
    if args.explain:
        print("\n" + explain_partner(args.explain, sub_edges, edges.attrs))
        return

    # CLI preview
    print("\nTop preview:")
    print(ranked.head(min(args.top_k, 15)).to_string(index=False))


if __name__ == "__main__":
    main()
