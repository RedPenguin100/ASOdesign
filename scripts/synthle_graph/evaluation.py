
import pandas as pd

# Open the genetic csv (assume the filename as in synthle.py)

# df_genetic = pd.read_csv("tp53_genetic_results.csv")
# df_physical = pd.read_csv("tp53_physical_results.csv")
# df_experimental = pd.read_csv("tp53_experimental_results.csv")
# df_db = pd.read_csv("tp53_database_results.csv")
# df_coexp = pd.read_csv("tp53_coexpression_results.csv")

# For each dataframe, rename columns (except 'CandidateGene') to include the source
# df_genetic = df_genetic.rename(columns={col: f"{col}_genetic" for col in df_genetic.columns if col != "CandidateGene"})
# df_physical = df_physical.rename(columns={col: f"{col}_physical" for col in df_physical.columns if col != "CandidateGene"})
# df_experimental = df_experimental.rename(columns={col: f"{col}_experimental" for col in df_experimental.columns if col != "CandidateGene"})
# df_db = df_db.rename(columns={col: f"{col}_db" for col in df_db.columns if col != "CandidateGene"})
# df_coexp = df_coexp.rename(columns={col: f"{col}_coexp" for col in df_coexp.columns if col != "CandidateGene"})

# df_genetic.to_csv("tp53_genetic_results.csv", index=False)
# df_physical.to_csv("tp53_physical_results.csv", index=False)
# df_experimental.to_csv("tp53_experimental_results.csv", index=False)
# df_db.to_csv("tp53_database_results.csv", index=False)
# df_coexp.to_csv("tp53_coexpression_results.csv", index=False)


# df_synth = df_synth.rename(columns={"second_row": "CandidateGene"})
# df_synth = df_synth.rename(columns={"third_row": "SLIDR mut pvalue"})
# df_synth = df_synth.rename(columns={"forth_row": "r.statistic_score"})

# # df_synth = df_synth.iloc[1:].reset_index(drop=True)
# df_synth.to_csv("synthetic_candidate.csv", index=False)




# df_synth.to_csv("synthetic_candidate.csv", index=False)
import numpy as np
df_merged = pd.read_csv("merged_updated.csv")
df_merged = df_merged.replace(r'^\s*$', np.nan, regex=True)

# Save the DataFrame as an Excel file, preserving NaN for unwritten cells
df_merged.to_excel("merged_updated.xlsx", index=False, na_rep="NaN")

# df_synth = pd.read_csv("synthetic_candidate.csv")

# # Perform inner join on 'CandidateGene', with df_synth columns first
# df_merged = df_synth.merge(df_merged, on="CandidateGene", how="inner")

# print(df_merged.head())
# df_merged.to_csv("merged_updated.csv", index=False)

# Outer join all dataframes on the CandidateGene column
# merged = df_genetic.merge(df_physical, on="CandidateGene", how="outer") \
#                .merge(df_experimental, on="CandidateGene", how="outer") \
#                .merge(df_db, on="CandidateGene", how="outer") \
#                .merge(df_coexp, on="CandidateGene", how="outer")

# Save the merged dataframe to a new CSV file
# merged.to_csv("merged.csv", index=False)





