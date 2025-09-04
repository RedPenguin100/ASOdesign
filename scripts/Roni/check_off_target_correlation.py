import pandas as pd
import matplotlib.pyplot as plt
import glob
import re

# --- SETTINGS ---
may_df = pd.read_csv("/home/oni/ASOdesign/scripts/data_genertion/data_updated_inhibition.csv")
inhibition = may_df[['index', 'Inhibition(%)']]

# file naming scheme: off_target.GEO.top{N}.cutoff{C}.premRNA.csv
pattern = "/home/oni/ASOdesign/scripts/data_genertion/off_target/off_target.GEO.top*.cutoff*.premRNA.csv"

# color mapping
color_map = {
    "50": "orange",
    "75": "purple",
    "100": "blue",
    "150": "cyan"
}

plt.figure(figsize=(10, 6))

for file in glob.glob(pattern):
    # extract topN and cutoff from filename
    m = re.search(r'top(\d+)\.cutoff(\d+)', file)
    if not m:
        continue
    topN, cutoff = m.groups()

    df = pd.read_csv(file)
    df = df[['index', 'off_target_score_MS']]
    merged = pd.merge(inhibition, df, on="index", how="inner")

    # correlation including zeros
    corr_all = merged['Inhibition(%)'].corr(merged['off_target_score_MS'])
    plt.scatter(int(cutoff), corr_all,
                color=color_map.get(topN, "gray"),
                marker='o',
                label=f"top {topN}" if f"top {topN}" not in plt.gca().get_legend_handles_labels()[1] else "")

    # correlation excluding zeros
    merged_nz = merged[merged['off_target_score_MS'] != 0]
    if len(merged_nz) > 1:
        corr_nz = merged_nz['Inhibition(%)'].corr(merged_nz['off_target_score_MS'])
        plt.scatter(int(cutoff), corr_nz,
                    color=color_map.get(topN, "gray"),
                    marker='x',
                    label=f"top {topN}" if f"top {topN}" not in plt.gca().get_legend_handles_labels()[1] else "")

plt.xlabel("Cutoff")
plt.ylabel("Correlation with Inhibition(%)")
plt.title("Correlation of Off-target Scores with Inhibition")
plt.legend(title="topN", loc="best")
plt.grid(True)
plt.show()