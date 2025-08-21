import pandas as pd
import matplotlib.pyplot as plt
import glob
import re

# --- SETTINGS ---
may_df = pd.read_csv("/home/oni/ASOdesign/scripts/data_genertion/data_updated_inhibition.csv")
inhibition = may_df[['index', 'Inhibition(%)']]

# file naming scheme: off_target.top{N}.cutoff{C}.premRNA.csv
pattern = "/home/oni/ASOdesign/scripts/data_genertion/off_target.GEO.top*.cutoff*.premRNA.csv"

# assign a color per "topN"
color_map = {
    "50": "blue",
    "75": "green",
    "100": "orange",
    "150": "red"
}
# shades per score variant
shade_styles = {
    "MS": "o",           # original log
    "MS_nozero": "x",    # log with zeros removed
}

results = []

for filepath in glob.glob(pattern):
    # extract numbers from filename
    m = re.search(r"top(\d+)\.cutoff(\d+)", filepath)
    if not m:
        continue
    topN, cutoff = m.groups()

    off_target_df = pd.read_csv(filepath)

    for col in ["off_target_score_MS"]:
        merged = pd.merge(inhibition, off_target_df[['index', col]], on="index", how="inner")

        # normal correlation
        corr = merged["Inhibition(%)"].corr(merged[col])
        results.append({
            "topN": topN,
            "cutoff": int(cutoff),
            "variant": col.split("_")[-1],  # log / TPM
            "corr": corr
        })

        # correlation after removing zeros
        filtered = merged[merged[col] != 0]
        if len(filtered) > 1:
            corr_nozero = filtered["Inhibition(%)"].corr(filtered[col])
            results.append({
                "topN": topN,
                "cutoff": int(cutoff),
                "variant": col.split("_")[-1] + "_nozero",
                "corr": corr_nozero
            })

# put results in dataframe
res_df = pd.DataFrame(results)

# --- PLOTTING ---
plt.figure(figsize=(10, 6))

for topN, group in res_df.groupby("topN"):
    base_color = color_map[topN]
    for variant, subg in group.groupby("variant"):
        plt.plot(
            subg["cutoff"], subg["corr"],
            marker=shade_styles[variant],
            linestyle="",
            label=f"Top {topN} - {variant}" if variant == "log" else None,  # avoid duplicate legend spam
            color=base_color,
            alpha=0.6 if "nozero" in variant else 1.0
        )

plt.xlabel("Cutoff")
plt.ylabel("Correlation (Inhibition vs Off-target score)")
plt.title("{Off-target X Inhibition} correlation across cutoffs and topN genes")
plt.grid(True)

# Custom legend (one entry per topN, plus markers explained separately)
from matplotlib.lines import Line2D
legend_elements = []
for topN, color in color_map.items():
    legend_elements.append(Line2D([0], [0], color=color, marker='o', linestyle="", label=f"Top {topN}"))

for variant, marker in shade_styles.items():
    legend_elements.append(Line2D([0], [0], color="black", marker=marker, linestyle="", label=variant))

plt.legend(handles=legend_elements, title="Legend", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()