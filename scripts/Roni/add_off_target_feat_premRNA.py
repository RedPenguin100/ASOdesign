import pandas as pd
import os

# =============================================== Off-target Function =================================================
from add_off_target_feat import get_off_target_feature



# =============================================== Load main ASO data ==================================================
script_dir = os.path.dirname(__file__)
data_dir = os.path.abspath(os.path.join(script_dir, "..", "data_genertion"))

data_path = os.path.join(data_dir, "data_updated_inhibition.csv")
may_df = pd.read_csv(data_path)
#may_df = may_df.head(1)

path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/'

# Using the mutated pre-mRNA sequences
sequences = '.mutated_transcriptome_premRNA.merged.csv'

cell_line_list = ['ACH-001328', 'ACH-000463', 'ACH-001188', 'ACH-001086', 'ACH-000739', 'ACH-000232']

cellline_to_filename = {
    'A431': 'ACH-001328',
    'A-431': 'ACH-001328',
    'NCI-H460': 'ACH-000463',
    'SH_SY5Y': 'ACH-001188',
    'HeLa': 'ACH-001086',
    'Hela': 'ACH-001086',
    'HepG2': 'ACH-000739',
    'U-251MG': 'ACH-000232'
}

# Final storage for all scores
all_results_TPM = {}
all_results_norm = {}


for cell_line, group_df in may_df.groupby('Cell_line'):
    print(f'Processing {cell_line}...')

    # Process only for cell lines which we have their sequences
    if cell_line not in cellline_to_filename:
        print(f"Skipping {cell_line}, no file mapping.")
        continue

    # Loading the relevant Sequences Data
    file_id = cellline_to_filename[cell_line]
    transcriptome_path = path+file_id+sequences
    if not os.path.exists(transcriptome_path):
        print(f"File missing for {cell_line}: {transcriptome_path}")
        continue

    print(f"Processing cell line {cell_line}")
    curr_df = pd.read_csv(transcriptome_path)
    print(f"Read {transcriptome_path}")

    # cut in case too heavy
    n = 500
    curr_df = curr_df.head(n)

    cell_line2df_ = {cell_line: curr_df}

    scores_TPM, scores_norm = get_off_target_feature(cell_line2df_, group_df)
    print(f"Done processing {cell_line}")

    all_results_TPM.update(scores_TPM)
    all_results_norm.update(scores_norm)

# Combine into final DataFrame
off_target_feature = pd.DataFrame({
    "index": list(all_results_TPM.keys()),
    "off_target_score_TPM": list(all_results_TPM.values()),
    "off_target_score_log": list(all_results_norm.values())
})

out_path = os.path.join(data_dir, "off_target_feature_500_premRNA.csv")
off_target_feature.to_csv(out_path, index=False)
print(f"Saved {out_path}")