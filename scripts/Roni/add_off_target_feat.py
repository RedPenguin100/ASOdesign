from asodesigner.fold import get_trigger_mfe_scores_by_risearch
from scripts.Roni.off_target_functions import dna_to_rna_reverse_complement, parse_risearch_output, aggregate_off_targets
from mutate_cell_line_transcriptome import celline_list

import pandas as pd
from io import StringIO
import os

'''

path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/'
sequences = '.mutated_transcriptome_premRNA.merged.csv'



# Set base directory relative to script
script_dir = os.path.dirname(__file__)
data_dir = os.path.abspath(os.path.join(script_dir, "..", "data_genertion"))

# ========================================== Pre-Processing Main ASO data =============================================
# Load the main ASO dataset
data_path = os.path.join(data_dir, "data_updated_inhibition.csv")
may_df = pd.read_csv(data_path)
may_df = may_df.head(1)

# Expression files path

# Load each transcriptome file
A431_df = pd.read_csv(path+celline_list[0]+sequences)
print("1")
NCI_H460_df = pd.read_csv(path+celline_list[1]+sequences)
print("2")
SH_SY5Y_df = pd.read_csv(path+celline_list[2]+sequences)
print("3")
HeLa_df = pd.read_csv(path+celline_list[3]+sequences)
print("4")
HepG2_df = pd.read_csv(path+celline_list[4]+sequences)
print("5")
U_251MG_df = pd.read_csv(path+celline_list[5]+sequences)
print("6")

# ============================ Cut to top n expressed ==============================
n = 1
A431_df = A431_df.head(n)
NCI_H460_df = NCI_H460_df.head(n)
SH_SY5Y_df = SH_SY5Y_df.head(n)
HeLa_df = HeLa_df.head(n)
HepG2_df = HepG2_df.head(n)
U_251MG_df = U_251MG_df.head(n)
 # =============================================== Cell lines  ========================================================
cell_line_list = ['ACH-001328', 'ACH-000463', 'ACH-001188', 'ACH-001086', 'ACH-000739', 'ACH-000232']

cell_line2df_ = {'A431':A431_df,
                'A-431':A431_df,
                'NCI-H460':NCI_H460_df,
                'SH_SY5Y':SH_SY5Y_df,
                'HeLa':HeLa_df,
                'HepG2':HepG2_df,
                'U-251MG':U_251MG_df}
'''
# =============================================== Function ============================================================


def get_off_target_feature(cell_line2df, ASO_df):
    index_score_vec_TPM = {}
    index_score_vec_norm = {}

    for idx, row in ASO_df.iterrows():

        index = row['index']
        cell_line = row['Cell_line']
        ASO_seq = row['Sequence']
        trigger = dna_to_rna_reverse_complement(ASO_seq)
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
            mut_seq = gene_row.get('Mutated Transcript Sequence')
            og_seq = gene_row.get('Original Transcript Sequence')

            if pd.isna(mut_seq) and pd.isna(og_seq):
                continue

            mRNA_seq = mut_seq if not pd.isna(mut_seq) else og_seq
            exp_TPM = gene_row['expression_TPM']
            exp_norm = gene_row['expression_norm']

            name_to_seq[curr_gene] = mRNA_seq
            name_to_exp_TPM[curr_gene] = exp_TPM
            name_to_exp_norm[curr_gene] = exp_norm

        if not name_to_seq:
            print(f"[WARNING] name_to_seq is empty for index={index} cell_line={cell_line}")

        # Calculate mfe scores
        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
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

'''
off_target_vec = get_off_target_feature(cell_line2df_, may_df)
#print(off_target_vec)

off_target_feature = pd.DataFrame({
    "off_target_score_TPM": off_target_vec[0],
    "off_target_score_log": off_target_vec[1],
}).reset_index().rename(columns={"index": "index"})

off_target_feature.to_csv(os.path.join(data_dir, "off_target.top500.cutoff1200.premRNA.csv"), index=False)
print("off_target_feature.csv saved")

'''