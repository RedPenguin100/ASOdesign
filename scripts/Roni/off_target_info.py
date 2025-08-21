from ASOdesign.scripts.Roni.add_off_target_feat import NCI_H460_df, U_251MG_df
from ASOdesign.scripts.Roni.off_target_info import SH_SY5Y_df, HepG2_df
from asodesigner.fold import get_trigger_mfe_scores_by_risearch
from scripts.Roni.off_target_functions import dna_to_rna_reverse_complement, parse_risearch_output, aggregate_off_targets
import pandas as pd
from io import StringIO
from asodesigner.consts import DATA_PATH_NEW
import os
import pickle



# Set base directory relative to script
script_dir = os.path.dirname(__file__)
data_dir = os.path.abspath(os.path.join(script_dir, "..", "data_genertion"))

# Load the main ASO dataset
data_path = os.path.join(data_dir, "data_asoptimizer_updated.csv")
may_df = pd.read_csv(data_path)
may_df = may_df.head(10)

# Expression files path
expr_path = os.path.join(data_dir, "cell_line_expression")

# Load each transcriptome file -================== MATURE mRNA ===========================
'''
A431_df = pd.read_csv(os.path.join(expr_path, 'ACH-001328_transcriptome.csv'))
NCI_H460_df = pd.read_csv(os.path.join(expr_path, 'ACH-000463_transcriptome.csv'))
SH_SY5Y_df = pd.read_csv(os.path.join(expr_path, 'ACH-001188_transcriptome.csv'))
HeLa_df = pd.read_csv(os.path.join(expr_path, 'ACH-001086_transcriptome.csv'))
HepG2_df = pd.read_csv(os.path.join(expr_path, 'ACH-000739_transcriptome.csv'))
U_251MG_df = pd.read_csv(os.path.join(expr_path, 'ACH-000232_transcriptome.csv'))
'''
# Load each transcriptome file -================== PRE mRNA ===========================
# A431_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-001328.mutated_transcriptome_premRNA.merged.csv')
# NCI_H460_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-000463.mutated_transcriptome_premRNA.merged.csv')
# SH_SY5Y_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-001188.mutated_transcriptome_premRNA.merged.csv')
# HeLa_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-001086.mutated_transcriptome_premRNA.merged.csv')
# HepG2_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-000739.mutated_transcriptome_premRNA.merged.csv')
# U_251MG_df = pd.read_csv(DATA_PATH_NEW/'cell_line_expression/ACH-000232.mutated_transcriptome_premRNA.merged.csv')





# ============================ Cut to top n expressed ==============================


# A431_df = A431_df.head(n)
# NCI_H460_df = NCI_H460_df.head(n)
# SH_SY5Y_df = SH_SY5Y_df.head(n)
# HeLa_df = HeLa_df.head(n)
# HepG2_df = HepG2_df.head(n)
# U_251MG_df = U_251MG_df.head(n)
 # ===================================================================================
cell_line_list = ['ACH-001328', 'ACH-000463', 'ACH-001188', 'ACH-001086', 'ACH-000739', 'ACH-000232']

# cell_line2df_ = {
#     'A431':A431_df,
#     'A-431':A431_df,
#     'NCI-H460':NCI_H460_df,
#     'SH_SY5Y':SH_SY5Y_df,
#     'SH-SY5Y':SH_SY5Y_df,
#     'HeLa':HeLa_df,
#     'Hela':HeLa_df,
#     'HepG2':HepG2_df,
#     'U-251MG':U_251MG_df,
#     'U251':U_251MG_df
# }

cellline_to_filename = {
    'A431': 'ACH-001328',
    'A-431': 'ACH-001328',
    'NCI-H460': 'ACH-000463',
    'SH_SY5Y': 'ACH-001188',
    'SH-SY5Y': 'ACH-001188',
    'HeLa': 'ACH-001086',
    'Hela': 'ACH-001086',
    'HepG2': 'ACH-000739',
    'U-251MG': 'ACH-000232',
    'U251': 'ACH-000232'
}


def get_off_target_info(cell_line2df, ASO_df, cutoff, topN):
    index_info_vec = {}

    for idx, row in ASO_df.iterrows():

        index = row['index']
        cell_line = row['Cell_line']
        ASO_seq = row['Sequence']
        trigger = dna_to_rna_reverse_complement(ASO_seq)
        target_gene = row['Canonical Gene Name']

        if cell_line not in cell_line2df:
            continue

        cell_line_id = cell_line2df[cell_line]
        curr_df = pd.read_csv(DATA_PATH_NEW/f'cell_line_expression/{cell_line_id}.mutated_transcriptome_premRNA.merged.csv')
        curr_df = curr_df.head(topN)
        print(f'read {cell_line_id} data')

        name_to_seq = {}
        name_to_exp_TPM = {}
        name_to_exp_norm = {}
        #name_to_transcript_id = {}

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
            transcript_id = gene_row['Transcript_ID']

            name_to_seq[curr_gene] = mRNA_seq
            name_to_exp_TPM[curr_gene] = exp_TPM
            name_to_exp_norm[curr_gene] = exp_norm
            #name_to_transcript_id[curr_gene] = transcript_id

        # Calculate mfe scores
        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
            name_to_seq,
            minimum_score=cutoff,
            parsing_type='2'
        )

        #print(f'000:\n{result_dict}\n')

        result_df = parse_risearch_output(result_dict)

        top_score_df = result_df.loc[result_df.groupby('target')['score'].idxmax()]
        top_score_df = top_score_df.sort_values(by='score', ascending=False).reset_index(drop=True)
        #print(f'001:\n{result_df}\n{result_df.columns}')
        #print(f'001:\n{top_score_df}\n{top_score_df.columns}')

        top_score_df['energy_weighted_TPM'] = top_score_df['energy'] * top_score_df['target'].map(name_to_exp_TPM)
        top_score_df['energy_weighted_norm'] = top_score_df['energy'] * top_score_df['target'].map(name_to_exp_norm)

        #print(f'002:\n{top_score_df}\n{top_score_df.columns}')
        index_info_vec[index] = top_score_df


    return index_info_vec

n = 50
cutoff = 1200

off_target_info_vec = get_off_target_info(cellline_to_filename, may_df, cutoff, n)
with open(f'off_target_info.premRNA.top{n}.cutoff{cutoff}TEST.pkl', 'wb') as f:
    pickle.dump(off_target_info_vec, f)
print("yay saved")
# off_target_feature = pd.DataFrame({
#     "off_target_score_TPM": off_target_vec[0],
# }).reset_index().rename(columns={"index": "index"})

#off_target_feature.to_csv(os.path.join(data_dir, "off_target___.csv"), index=False)
#print("off_target_feature.csv saved")

