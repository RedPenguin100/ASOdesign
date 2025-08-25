###### imports ######
from typing import Dict, Optional, Tuple, List
import pandas as pd
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt

###### helpers ######
def normalize_enst(x):
    """
    Remove ENST version suffix.
    Example: 'ENST00000331789.3' -> 'ENST00000331789'
    input:
        x: any value convertible to string.
    Returns:
        str: normalized ENST without the version suffix.
    """
    s = str(x).strip()
    if "." in s:
        return s.split(".")[0]
    return s


def gene_symbol_from_gene_col(g):
    """
    Extract the gene symbol from a 'Gene' column value.
    Example: 'TMSB4X (7114)' -> 'TMSB4X'
    input:
        g: any value convertible to string (e.g., 'TMSB4X (7114)').
    Returns:
        str: the gene symbol without trailing '(...)'.
    """
    s = str(g).strip()
    if "(" in s:
        return s.split("(")[0].strip()
    return s


def find_premrna_for_target(cell_df, target):
    """
    Find the pre-mRNA sequence for a given target (ENST or gene symbol)
    inside a cell-line dataframe.

    Expects 'cell_df' to have at least these columns:
      - 'Transcript_ID'
      - 'Gene'
    And (preferably) one of the sequence columns:
      - 'Mutated Transcript Sequence'
      - 'Original Transcript Sequence'

    inputs:
        cell_df: pandas DataFrame of the cell line transcriptome.
        target:  ENST ID (with or without version) or a gene symbol (e.g., 'KRAS').

    Returns:
        (sequence_string, source_column_name) if found,
        otherwise (None, None).
    """
    if "Transcript_ID" not in cell_df.columns or "Gene" not in cell_df.columns:
        raise KeyError("Cell-line dataframe must have 'Transcript_ID' and 'Gene' columns.")

    # 1) Try match by Transcript_ID (normalized without version)
    enst_norm = normalize_enst(target)
    hits = cell_df[cell_df["Transcript_ID"].astype(str).map(normalize_enst) == enst_norm]

    # 2) Fallback: match by gene symbol (strip '(Entrez)' part)
    if len(hits) == 0:
        mask = cell_df["Gene"].astype(str).map(gene_symbol_from_gene_col) == str(target)
        hits = cell_df[mask]

    if len(hits) == 0:
        return None, None

    row = hits.iloc[0]

    # Prefer mutated sequence; fallback to original
    seq = row.get("Mutated Transcript Sequence")
    if isinstance(seq, str) and len(seq) > 0:
        return seq, "Mutated Transcript Sequence"

    seq = row.get("Original Transcript Sequence")
    if isinstance(seq, str) and len(seq) > 0:
        return seq, "Original Transcript Sequence"

    return None, None

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Try to import p-values from scipy, fall back gracefully if not installed
try:
    from scipy.stats import pearsonr, spearmanr
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

##### plot and analyse #####

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Try to import p-values from scipy, fall back gracefully if not installed
try:
    from scipy.stats import pearsonr, spearmanr
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


def corr_and_scatter_df(
    df_feature: pd.DataFrame,
    df_all_subset: pd.DataFrame,
    feature_col: str = None,                 # if None and df_feature has a single non-index col → auto-detect
    index_col: str = "index",
    inhib_col: str = "log_inhibition",
    remove_outliers: str = "iqr",            # "iqr", "percentile", or None
    iqr_k: float = 1.5,
    p_low: float = 0.01, p_high: float = 0.99,
    aggregate_feature: str = "max",          # "max" / "mean" / "median"
    aggregate_inhib: str = "median",         # "median" / "mean" / "max"
    title: str = None,
    figsize=(7, 6),
    save_path: str | None = None,            # if provided, saves the figure
):
    """
    Merge a feature DataFrame with labels (log_inhibition), compute Pearson/Spearman,
    optionally remove outliers, and plot a scatter with a simple regression line.

    Parameters
    ----------
    df_feature : pd.DataFrame
        Must contain [index_col, <feature_col>].
    df_all_subset : pd.DataFrame
        Subset of df_all for a single cell line. Must contain [index_col, inhib_col].
    feature_col : str | None
        Name of the feature column. If None and df_feature has exactly one non-index column,
        it will be auto-detected.
    """

    # --- basic checks / copy ---
    df_feat = df_feature.copy()
    if index_col not in df_feat.columns:
        raise KeyError(f"Feature table must contain '{index_col}'")

    # auto-detect feature column if needed
    if feature_col is None:
        non_idx_cols = [c for c in df_feat.columns if c != index_col]
        if len(non_idx_cols) != 1:
            raise ValueError("feature_col is None but could not auto-detect uniquely. "
                             f"Candidates: {non_idx_cols}")
        feature_col = non_idx_cols[0]

    # --- aggregate duplicates in feature by index ---
    if aggregate_feature not in {"max","mean","median"}:
        raise ValueError("aggregate_feature must be one of {'max','mean','median'}")
    agg_fun_feat = {"max": "max", "mean": "mean", "median": "median"}[aggregate_feature]

    df_feat_agg = (df_feat[[index_col, feature_col]]
                   .dropna(subset=[feature_col])
                   .groupby(index_col, as_index=False)
                   .agg({feature_col: agg_fun_feat}))

    # --- aggregate duplicates in labels by index ---
    if index_col not in df_all_subset.columns or inhib_col not in df_all_subset.columns:
        raise KeyError(f"df_all_subset must contain '{index_col}' and '{inhib_col}'")

    if aggregate_inhib not in {"median","mean","max"}:
        raise ValueError("aggregate_inhib must be one of {'median','mean','max'}")
    agg_fun_lab = {"median": "median", "mean": "mean", "max": "max"}[aggregate_inhib]

    labels = (df_all_subset[[index_col, inhib_col]]
              .dropna(subset=[inhib_col])
              .groupby(index_col, as_index=False)
              .agg({inhib_col: agg_fun_lab}))

    # --- merge feature & labels ---
    merged = df_feat_agg.merge(labels, on=index_col, how="inner")

    # coerce to numeric
    x = pd.to_numeric(merged[feature_col], errors="coerce")
    y = pd.to_numeric(merged[inhib_col],   errors="coerce")
    mask = x.notna() & y.notna()
    merged = merged.loc[mask].copy()

    # --- outlier removal ---
    def _apply_outliers(df, cols):
        if remove_outliers is None:
            return df
        if remove_outliers == "iqr":
            keep = pd.Series(True, index=df.index)
            for c in cols:
                q1, q3 = df[c].quantile([0.25, 0.75])
                iqr = q3 - q1
                lo, hi = q1 - iqr_k * iqr, q3 + iqr_k * iqr
                keep &= df[c].between(lo, hi)
            return df[keep].copy()
        if remove_outliers == "percentile":
            keep = pd.Series(True, index=df.index)
            for c in cols:
                lo, hi = df[c].quantile([p_low, p_high])
                keep &= df[c].between(lo, hi)
            return df[keep].copy()
        raise ValueError("remove_outliers must be 'iqr', 'percentile', or None")

    merged = _apply_outliers(merged, [feature_col, inhib_col])

    # recompute x, y after filtering
    x = pd.to_numeric(merged[feature_col], errors="coerce")
    y = pd.to_numeric(merged[inhib_col],   errors="coerce")
    mask = x.notna() & y.notna()
    merged = merged.loc[mask].copy()
    x = x[mask]; y = y[mask]

    n = len(merged)
    if n < 2:
        print("[WARN] Not enough points to compute correlation.")
        stats = {"n": n, "pearson_r": None, "pearson_p": None,
                 "spearman_r": None, "spearman_p": None,
                 "feature_col": feature_col}
        return stats, merged

    # --- correlations ---
    if _HAS_SCIPY:
        pr, pp = pearsonr(x, y)
        sr, sp = spearmanr(x, y)
    else:
        pr, pp = x.corr(y, method="pearson"), None
        sr, sp = x.corr(y, method="spearman"), None

    print(f"n = {n}")
    print(f"Pearson r = {pr:.3f}" + (f", p = {pp:.3g}" if pp is not None else " (no p-value; SciPy not installed)"))
    print(f"Spearman r = {sr:.3f}" + (f", p = {sp:.3g}" if sp is not None else " (no p-value; SciPy not installed)"))

    # --- scatter plot ---
    plt.figure(figsize=figsize)
    plt.scatter(x, y, s=20, alpha=0.5)
    if n >= 2:
        try:
            m, b = np.polyfit(x, y, 1)
            xx = np.linspace(x.min(), x.max(), 200)
            yy = m * xx + b
            plt.plot(xx, yy)
        except Exception:
            pass

    ttl = title or f"Feature vs {inhib_col}"
    plt.title(ttl)
    plt.xlabel(feature_col)
    plt.ylabel(inhib_col)
    plt.grid(True, alpha=0.3)

    # annotate correlations
    txt = f"Pearson r={pr:.3f}"
    if pp is not None:
        txt += f", p={pp:.1g}"
    txt += f"\nSpearman r={sr:.3f}"
    if sp is not None:
        txt += f", p={sp:.1g}"
    plt.gca().text(0.02, 0.98, txt, transform=plt.gca().transAxes, va="top", ha="left")

    if save_path:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"[OK] saved figure → {os.path.abspath(save_path)}")
    else:
        plt.show()

    stats = {"n": int(n),
             "pearson_r": float(pr), "pearson_p": (float(pp) if pp is not None else None),
             "spearman_r": float(sr), "spearman_p": (float(sp) if sp is not None else None),
             "feature_col": feature_col}
    return stats, merged



####### main #######
def compute_offtarget_folding_feature(
    aso_pkl_path: str,
    df_all: pd.DataFrame,                       # must contain ASO id ('index') and a cell-line column
    cell_line2df: Dict[str, pd.DataFrame],      # cell line name -> transcriptome dataframe
    *,
    score_col: str = 'energy_w_by_exp_ARTH',    # which PKL column to use
    use_abs_for_selection: bool = True,         # choose off-target by max |score| if True
    use_abs_in_feature: bool = True,            # feature uses |score| if True (else signed)
    flank_size: int = 120,
    window_size: int = 45,
    step: int = 7,
    coords_are_1based_inclusive: bool = True,   # PKL coordinates start at 1 and are inclusive
    eps_min_abs_mfe: float = 1e-6,              # avoid division by zero
    output_csv: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    For each ASO:
      1) pick the strongest off-target row (by `score_col`)
      2) get its cell line from `df_all`
      3) in the matching cell-line table, resolve 'target' to a transcript row
      4) extract pre-mRNA sequence (mutated if available, else original)
      5) compute avg MFE over the sense region within a flanked context
      6) feature = (chosen score) / |MFE|
    Returns a tidy DataFrame; optionally writes CSV.
    """
    # --- pick the cell-line column from df_all ---
    if "cell_line_uniform" in df_all.columns:
        cell_col = "cell_line_uniform"
    elif "Cell_line" in df_all.columns:
        cell_col = "Cell_line"
    else:
        raise KeyError("df_all must contain either 'cell_line_uniform' or 'Cell_line'.")

    if "index" not in df_all.columns:
        raise KeyError("df_all must contain the ASO identifier column 'index'.")

    # --- load and flatten PKL: {aso_idx -> DataFrame} -> single DataFrame ---
    if not os.path.exists(aso_pkl_path):
        raise FileNotFoundError(f"PKL not found: {aso_pkl_path}")
    with open(aso_pkl_path, "rb") as f:
        store = pickle.load(f)

    frames: List[pd.DataFrame] = []
    for aso_idx, ds in store.items():
        if isinstance(ds, pd.DataFrame) and not ds.empty:
            tmp = ds.copy()
            tmp.insert(0, "aso_idx", int(aso_idx))
            frames.append(tmp)

    # If nothing to do, return/save empty dataframe
    if not frames:
        cols = ["aso_idx", "cell_line", "target", "target_start", "target_end",
                score_col, "mfe_offtarget", "feature", "seq_col"]
        out = pd.DataFrame(columns=cols)
        if output_csv:
            out.to_csv(output_csv, index=False)
            if verbose:
                print(f"Saved (empty): {os.path.abspath(output_csv)}")
        return out

    df_off = pd.concat(frames, ignore_index=True)

    # --- required columns in off-target table ---
    req_cols = {"aso_idx", "target", "target_start", "target_end", score_col}
    missing = req_cols - set(df_off.columns)
    if missing:
        raise KeyError(f"Off-target table missing columns: {missing}")

    # Ensure numeric score; drop rows with NaN scores
    df_off[score_col] = pd.to_numeric(df_off[score_col], errors="coerce")
    df_off = df_off.dropna(subset=[score_col])
    if df_off.empty:
        cols = ["aso_idx", "cell_line", "target", "target_start", "target_end",
                score_col, "mfe_offtarget", "feature", "seq_col"]
        out = pd.DataFrame(columns=cols)
        if output_csv:
            out.to_csv(output_csv, index=False)
            if verbose:
                print(f"Saved (empty after score filter): {os.path.abspath(output_csv)}")
        return out

    # --- choose the strongest off-target per ASO ---
    sel_metric = df_off[score_col].abs() if use_abs_for_selection else df_off[score_col]
    df_off["_sel_metric"] = sel_metric
    top_idx = df_off.groupby("aso_idx")["_sel_metric"].idxmax()
    top = df_off.loc[top_idx, ["aso_idx", "target", "target_start", "target_end", score_col]].copy()

    # --- attach cell line per ASO from df_all ---
    idx2cl = (
        df_all[["index", cell_col]]
        .dropna()
        .drop_duplicates()
        .rename(columns={"index": "aso_idx", cell_col: "cell_line"})
    )
    # be robust to index dtype
    idx2cl["aso_idx"] = idx2cl["aso_idx"].astype(int)
    top = top.merge(idx2cl, on="aso_idx", how="left")

    # --- compute MFE per ASO (one row per ASO in `top`) ---
    from Folding_Functions import get_sense_with_flanks, calculate_avg_mfe_over_sense_region  # ensure imported

    mfes: List[float] = []
    seq_cols_used: List[Optional[str]] = []
    skipped = {"no_cell_line": 0, "no_cl_df": 0, "no_target_seq": 0, "bad_bounds": 0}

    for _, r in top.iterrows():
        # a) resolve cell-line table
        cl = r["cell_line"]
        if pd.isna(cl):
            mfes.append(np.nan); seq_cols_used.append(None); skipped["no_cell_line"] += 1; continue
        cl_name = str(cl).strip()
        cl_df = cell_line2df.get(cl_name)
        if cl_df is None:
            mfes.append(np.nan); seq_cols_used.append(None); skipped["no_cl_df"] += 1; continue

        # b) find pre-mRNA sequence for this target in that table
        seq_str, seq_col = find_premrna_for_target(cl_df, r["target"])
        if seq_str is None:
            mfes.append(np.nan); seq_cols_used.append(None); skipped["no_target_seq"] += 1; continue

        # c) convert coordinates (1-based inclusive → 0-based end-exclusive)
        if coords_are_1based_inclusive:
            start = int(r["target_start"]) - 1
            end   = int(r["target_end"])
        else:
            start = int(r["target_start"])
            end   = int(r["target_end"])

        # d) normalize sequence alphabet to RNA (U instead of T)
        seq_str = str(seq_str).upper().replace("T", "U")

        # e) sanity check bounds
        if not (0 <= start < end <= len(seq_str)):
            mfes.append(np.nan); seq_cols_used.append(seq_col); skipped["bad_bounds"] += 1; continue

        # f) compute average MFE over the sense region within a flanked context
        sense_len = end - start
        flanked = get_sense_with_flanks(
            pre_mrna=seq_str,
            sense_start=start,
            sense_length=sense_len,
            flank_size=flank_size,
        )
        mfe = calculate_avg_mfe_over_sense_region(
            sequence=flanked,
            sense_start=start,
            sense_length=sense_len,
            flank_size=flank_size,
            window_size=window_size,
            step=step,
        )
        mfes.append(mfe)
        seq_cols_used.append(seq_col)

    # attach computed columns
    top["mfe_offtarget"] = mfes
    top["seq_col"] = seq_cols_used

    # --- final feature ---
    top["mfe_abs"] = top["mfe_offtarget"].abs().clip(lower=eps_min_abs_mfe)  # avoid /0
    chosen_score = top[score_col].abs() if use_abs_in_feature else top[score_col]
    feature_name = f"{score_col}_overAbsMFE" if use_abs_in_feature else f"{score_col}_signed_overAbsMFE"
    top[feature_name] = chosen_score / top["mfe_abs"]

    # tidy output
    out_cols = ["aso_idx", "cell_line", "target", "target_start", "target_end",
                score_col, "mfe_offtarget", feature_name, "seq_col"]
    out = top[out_cols].copy()
    out = out.rename(columns={"aso_idx": "index"})
    # save csv
    if output_csv:
        out.to_csv(output_csv, index=False)
        if verbose:
            print(f"Saved: {os.path.abspath(output_csv)}")

    if verbose:
        print(
            f"Computed feature for {out[feature_name].notna().sum()} ASOs | "
            f"skips: {skipped}"
        )

    return out