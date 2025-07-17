from asodesigner.fold import get_trigger_mfe_scores_by_risearch

import pandas as pd
from io import StringIO
import os



# Set base directory relative to script
script_dir = os.path.dirname(__file__)
data_dir = os.path.abspath(os.path.join(script_dir, "..", "data_genertion"))

# Load the main ASO dataset
data_path = os.path.join(data_dir, "data_asoptimizer_updated.csv")
may_df = pd.read_csv(data_path)
may_df = may_df.head(100)

# Expression files path
expr_path = os.path.join(data_dir, "cell_line_expression")

# Load each transcriptome file
A431_df = pd.read_csv(os.path.join(expr_path, 'ACH-001328_transcriptome.csv'))
NCI_H460_df = pd.read_csv(os.path.join(expr_path, 'ACH-000463_transcriptome.csv'))
SH_SY5Y_df = pd.read_csv(os.path.join(expr_path, 'ACH-001188_transcriptome.csv'))
HeLa_df = pd.read_csv(os.path.join(expr_path, 'ACH-001086_transcriptome.csv'))
HepG2_df = pd.read_csv(os.path.join(expr_path, 'ACH-000739_transcriptome.csv'))
U_251MG_df = pd.read_csv(os.path.join(expr_path, 'ACH-000232_transcriptome.csv'))

# ============================ Cut to top n expressed ==============================
n = 500
A431_df = A431_df.head(n)
NCI_H460_df = NCI_H460_df.head(n)
SH_SY5Y_df = SH_SY5Y_df.head(n)
HeLa_df = HeLa_df.head(n)
HepG2_df = HepG2_df.head(n)
U_251MG_df = U_251MG_df.head(n)
 # ===================================================================================
cell_line_list = ['ACH-001328', 'ACH-000463', 'ACH-001188', 'ACH-001086', 'ACH-000739', 'ACH-000232']

cell_line2df_ = {'A431':A431_df,
                'A-431':A431_df,
                'NCI-H460':NCI_H460_df,
                'SH_SY5Y':SH_SY5Y_df,
                'HeLa':HeLa_df,
                'HepG2':HepG2_df,
                'U-251MG':U_251MG_df}


def parse_risearch_output(output_str: str) -> pd.DataFrame:
    columns = ["trigger", "trigger_start", "trigger_end", "target", "target_start", "target_end", "score", "energy"]
    df = pd.read_csv(StringIO(output_str.strip()), sep="\t", header=None, names=columns)
    return df

def aggregate_off_targets(df: pd.DataFrame) -> pd.DataFrame:
    # Aggregate: sum score (if it's hybridization hits) and take minimum (strongest) energy
    grouped = df.groupby("target").agg({
        "score": "sum",
        "energy": "min"
    }).reset_index()
    return grouped


def get_off_target_feature(cell_line2df, ASO_df):
    index_score_vec_TPM = {}
    index_score_vec_norm = {}

    for idx, row in ASO_df.iterrows():

        index = row['index']
        cell_line = row['Cell_line']
        ASO_seq = row['Sequence']
        target_gene = row['Canonical Gene Name']

        if cell_line not in cell_line2df:
            continue

        curr_df = cell_line2df[cell_line]

        name_to_seq = {}
        name_to_exp_TPM = {}
        name_to_exp_norm = {}

        for _, gene_row in curr_df.iterrows():
            curr_gene = gene_row['Gene'].split()[0]

            # Skip the target gene
            if curr_gene == target_gene:
                continue

            # Select mutated or original sequence
            mut_seq = gene_row.get('Mutated Transcript sequence')
            og_seq = gene_row.get('Original Transcript sequence')

            if pd.isna(mut_seq) and pd.isna(og_seq):
                continue

            mRNA_seq = mut_seq if not pd.isna(mut_seq) else og_seq
            exp_TPM = gene_row['expression_TPM']
            exp_norm = gene_row['expression_norm']

            name_to_seq[curr_gene] = mRNA_seq
            name_to_exp_TPM[curr_gene] = exp_TPM
            name_to_exp_norm[curr_gene] = exp_norm


        # Calculate mfe scores
        result_dict = get_trigger_mfe_scores_by_risearch(
            ASO_seq,
            name_to_seq,
            minimum_score=800,
            parsing_type='2'
        )
        result_df = parse_risearch_output(result_dict)
        #print(f'001:\n {result_df}')
        result_df_agg = aggregate_off_targets(result_df)
        #print(f'002:\n {result_df_agg}')

        # Weight by expression
        result_dict_weighted_TPM = {
            row['target']: row['energy'] * name_to_exp_TPM.get(row['target'], 0)
            for _, row in result_df_agg.iterrows()
        }
        result_dict_weighted_norm = {
            row['target']: row['energy'] * name_to_exp_norm.get(row['target'], 0)
            for _, row in result_df_agg.iterrows()
        }
        #print(result_dict_weighted_TPM)
        index_score_vec_TPM[index] = sum(result_dict_weighted_TPM.values())

        #print(result_dict_weighted_norm)
        index_score_vec_norm[index] = sum(result_dict_weighted_norm.values())


    return [index_score_vec_TPM, index_score_vec_norm]


off_target_vec = get_off_target_feature(cell_line2df_, may_df)
#print(off_target_vec)

off_target_feature = pd.DataFrame({
    "off_target_score_TPM": off_target_vec[0],
    "off_target_score_log": off_target_vec[1],
}).reset_index().rename(columns={"index": "index"})

off_target_feature.to_csv(os.path.join(data_dir, "off_target_feature_500.csv"), index=False)
print("off_target_feature.csv saved")

