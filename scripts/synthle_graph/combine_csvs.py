import pandas as pd
from pathlib import Path

# Load all rows
base = pd.read_csv("base.csv")
gen_high = pd.read_csv("tp53_like_results_genetic_high.csv")
gen_low  = pd.read_csv("tp53_like_results_genetic_low.csv")
phy_high = pd.read_csv("tp53_like_results_physical_high.csv")
phy_low  = pd.read_csv("tp53_like_results_physical_low.csv")

other = pd.read_csv('tp53_like_all.csv')

BASE_KEY = "CandidateGene"
OTHER_KEY = "Protein"

EDGE_KEY = "Gene"
TYPE_COL = "edge_type"

def merge_edge_type(df_all, edges_df):
    result = df_all.copy()
    if TYPE_COL not in edges_df.columns:
        # If no edge_type column, treat as a generic type
        subsets = [("generic", edges_df.copy())]
    else:
        subsets = []
        for etype, sub in edges_df.groupby(TYPE_COL):
            subsets.append((str(etype), sub.copy()))

    for etype, sub in subsets:
        # Rename non-key, non-type columns with suffix
        cols_to_rename = [c for c in sub.columns if c not in (EDGE_KEY, TYPE_COL)]
        rename_map = {c: f"{c}__{etype}" for c in cols_to_rename}
        sub_ren = sub.rename(columns=rename_map)

        # Merge on CandidateGene (left) == Gene (right)
        result = result.merge(sub_ren, left_on=BASE_KEY, right_on=EDGE_KEY, how="left")
        # Drop the duplicate key column from right
        result = result.drop(columns=[EDGE_KEY, TYPE_COL], errors="ignore")
    return result


def merge_with_label(base_df, other_df, label):
    # rename all non-key columns with suffix
    other_df = other_df.copy()
    rename_map = {c: f"{c}__{label}" for c in other_df.columns if c != OTHER_KEY}
    other_df = other_df.rename(columns=rename_map)
    merged = base_df.merge(other_df, left_on=BASE_KEY, right_on=OTHER_KEY, how="left")
    # drop the duplicated key column from the right
    merged = merged.drop(columns=[OTHER_KEY], errors="ignore")
    return merged


# Build outputs
genetic = merge_with_label(base, gen_high, "genetic_high")
genetic = merge_with_label(genetic, gen_low, "genetic_low")

physical = merge_with_label(base, phy_high, "physical_high")
physical = merge_with_label(physical, phy_low, "physical_low")

all_merged = merge_with_label(base, gen_high, "genetic_high")
all_merged = merge_with_label(all_merged, gen_low, "genetic_low")
all_merged = merge_with_label(all_merged, phy_high, "physical_high")
all_merged = merge_with_label(all_merged, phy_low, "physical_low")

all_merged = merge_edge_type(all_merged, other)


def format_sci(v):
    if pd.isna(v):
        return pd.NA
    try:
        f = float(v)
        return f"{f:.3e}"
    except Exception:
        return pd.NA

for col in all_merged.columns:
    if col not in ['CandidateGene']:
        all_merged[col] = pd.to_numeric(all_merged[col], errors="coerce").map(format_sci)

# Save
out_dir = Path(".")
paths = {
    "genetic": out_dir / "merged_base_plus_genetic.csv",
    "physical": out_dir / "merged_base_plus_physical.csv",
    "all": out_dir / "merged_base_plus_all.csv",
}
# genetic.to_csv(paths["genetic"], index=False)
# physical.to_csv(paths["physical"], index=False)
all_merged.to_csv(paths["all"], index=False)