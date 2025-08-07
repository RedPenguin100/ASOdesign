import pandas as pd
from alternative_splicing import *
from consts import *
from Bio.Seq import Seq
from asodesigner.consts import DATA_PATH_NEW



def add_snp_for_myh7():
    df_seq = pd.read_csv('Genes with Sequences.csv')
    mhy = df_seq[df_seq['gene'] == 'MYH7']['gene_sequence'].values[0]
    seq = 'CTGGGCTGGATGAGATCATTGCCAAGCTGACCAAGGAGA'
    location_org = mhy.find(seq)
    if location_org != -1:
        print('original seq was found in ', str(location_org))
        print('swapping for SNP in pre-mRNA in location ', str(location_org+19))
        seq_c = 'CTGGGCTGGATGAGATCATCGCCAAGCTGACCAAGGAGA'
        s = list(mhy)
        s[11981] = 'C'
        mhy_c = ''.join(s)
        location_new = mhy_c.find(seq_c)
        if location_new != -1:
            print('new seq with SNP was found in ', str(location_new))
            print('swap was made - T changed into C')
            new_row = pd.DataFrame([{'gene':'MYH7_SNP_715C','gene_id': None , 'gene_sequence': mhy_c}])
            df_seq = df_seq[df_seq['gene'] != 'MHY7_SNP_715C']
            df_seq_new = pd.concat([df_seq, new_row])
            df_seq_new.to_csv('Genes with Sequences.csv', index=False)
            print('updated file for gene with sequences')



def find_subseq_locations_with_mismatch_info(genes_df, asos_df, max_mismatches=1):
    """
    Searches each ASO sequence in its corresponding gene sequence allowing up to max_mismatches.
    Returns detailed info including the gene window, mismatch positions, and a reconciled sequence.

    Args:
        genes_df (pd.DataFrame): columns ['gene', 'gene_sequence']
        asos_df (pd.DataFrame): columns ['Sequence', 'gene name']
        max_mismatches (int): maximum allowed mismatches

    Returns:
        pd.DataFrame with columns:
            ['gene', 'Sequence', 'match_position', 'mismatches',
             'gene_window', 'aso_sequence', 'reconciled_from_aso', 'mismatch_positions']
    """
    results = []

    # Make sure gene names are strings
    #genes_df['gene'] = genes_df['gene'].astype(str)
    #aso_df[CANONICAL_GENE] = aso_df[CANONICAL_GENE].astype(str)

    # Convert gene data to dictionary for quick lookup
    gene_seq_dict = dict(zip(genes_df['gene'], genes_df['gene_sequence']))

    for _, row in asos_df.iterrows():
        aso = row[SEQUENCE]
        aso_seq = str(Seq(aso).reverse_complement())
        gene_name = row[CANONICAL_GENE]
        gene_seq = gene_seq_dict.get(gene_name)
        if gene_seq is None:
            continue  # no matching gene available

        sub_len = len(aso_seq)
        for i in range(len(gene_seq) - sub_len + 1):
            gene_window = gene_seq[i : i + sub_len]
            # find mismatches
            mismatch_positions = [idx for idx, (g, a) in enumerate(zip(gene_window, aso_seq)) if g != a]
            num_mismatches = len(mismatch_positions)
            if num_mismatches <= max_mismatches:
                # reconciled_from_aso: start from aso_seq, replace mismatched chars with gene's
                reconciled = list(aso_seq)
                for pos in mismatch_positions:
                    reconciled[pos] = gene_window[pos]
                reconciled_from_aso = ''.join(reconciled)  # should equal gene_window

                results.append({
                    'gene': gene_name,
                    'Sequence': aso,
                    'Reverse ASO': aso_seq,
                    'match_position': i,
                    'mismatches': num_mismatches,
                    'gene_window': gene_window,
                    'reconciled_from_aso': reconciled_from_aso,
                    'mismatch_positions': mismatch_positions,
                })

    return pd.DataFrame(results)

def check_function(df, df_with_info):
    # --- 1. Check Row Counts (most important initial check) ---
    print(f"\nOriginal df rows: {len(df)}")
    print(f"df_with_info rows after merge: {len(df_with_info)}")

    if len(df) == len(df_with_info):
        print("Row counts match: The merge did NOT duplicate rows from the original df.")
    else:
        print("WARNING: Row counts DO NOT match. The merge introduced/removed rows.")

    # --- 2. Check if 'index' is still unique in df_with_info ---
    print(f"\nIs 'index' unique in df (original)? {df['index'].is_unique}")
    print(f"Is 'index' unique in df_with_info (after merge)? {df_with_info['index'].is_unique}")

    if df_with_info['index'].is_unique:
        print("The 'index' remains unique in df_with_info, confirming no row duplication.")
    else:
        print("WARNING: 'index' is NOT unique in df_with_info. This indicates row duplication.")

    # --- 3. Verify 'index' values are preserved ---
    # You can check if all index values are still present and in the same order
    # (if the merge didn't reorder things, which a left merge usually tries to avoid).
    # This is typically true if the 'Sequence' column in df_with_info had unique values too.
    if df['index'].equals(df_with_info['index']):
        print("\n'index' values and order are identical between original df and df_with_info.")
    else:
        print(
            "\n'index' values/order might have changed. This is less critical if its uniqueness is maintained.")

    # --- 4. Inspect rows where new information was added (or not added) ---

    # Rows where new location information was successfully added
    print("\n--- Sample of rows where new location info was added (first 5) ---")
    # Filter for rows where Location_in_sequence is NOT NaN
    rows_with_location = df_with_info[df_with_info['Location_in_sequence'].notna()]
    print(rows_with_location.head())
    print(f"Total rows with location information: {len(rows_with_location)}")

    # Rows where no location information was found (i.e., new columns are NaN)
    print("\n--- Sample of rows where NO new location info was found (first 5) ---")
    # Filter for rows where Location_in_sequence IS NaN
    rows_without_location = df_with_info[df_with_info['Location_in_sequence'].isna()]
    print(rows_without_location.head())
    print(f"Total rows without location information: {len(rows_without_location)}")

    # --- 5. Compare specific values for a chosen 'Sequence' ---
    # Pick a Sequence that you know should have been mapped
    # Choose one from df_with_info that has a location, e.g., 'AGAGGAAAATTTTTCCATCAG'
    # And one that doesn't, if any.

    sample_sequence_with_loc = rows_with_location['Sequence'].iloc[0] if not rows_with_location.empty else None
    if sample_sequence_with_loc:
        print(f"\n--- Details for a sample sequence with location: '{sample_sequence_with_loc}' ---")
        print("From df (original):")
        print(df[df['Sequence'] == sample_sequence_with_loc])
        print("\nFrom df_with_info (after merge):")
        print(df_with_info[df_with_info['Sequence'] == sample_sequence_with_loc])

    # If you have a sequence that did NOT get a location, pick one:
    sample_sequence_without_loc = rows_without_location['Sequence'].iloc[0] if not rows_without_location.empty else None
    if sample_sequence_without_loc:
        print(f"\n--- Details for a sample sequence WITHOUT location: '{sample_sequence_without_loc}' ---")
        print("From df (original):")
        print(df[df['Sequence'] == sample_sequence_without_loc])
        print("\nFrom df_with_info (after merge):")
        print(df_with_info[df_with_info['Sequence'] == sample_sequence_without_loc])

    # --- 6. Check columns ---
    print("\nColumns in original df:", df.columns.tolist())
    print("Columns in df_with_info:", df_with_info.columns.tolist())
    print("New columns added:", [col for col in df_with_info.columns if col not in df.columns])



def find_and_mutate(sequences, subseq, mutate_pos_in_subseq=0):
    mutated_results  = []
    for seq_idx, seq in enumerate(sequences):
        new_seq = str(seq)  # start with original
        i = 0
        mutated = False
        while i <= len(new_seq) - len(subseq):
            if new_seq[i:i+len(subseq)] == subseq:
                target_index = i + mutate_pos_in_subseq
                if target_index < len(new_seq):
                    mutated_base = 'C'
                    new_seq = new_seq[:target_index] + mutated_base + new_seq[target_index+1:]
                    mutated = True
                i += len(subseq)
            else:
                i += 1
        if mutated:
            mutated_results.append((seq_idx, new_seq))

    return mutated_results


def change_gene_to_data_for_myh7(my_gene):
    genes_path = DATA_PATH_NEW / 'Genes with Sequences.csv'
    df_seq = pd.read_csv(genes_path)
    myh7_snp = df_seq[df_seq['gene'] == 'MYH7_SNP_715C']['gene_sequence'].values[0]
    subsequence = 'CTGGGCTGGATGAGATCATTGCCAAGCTGACCAAGGAGA'
    mutated_results_exons = find_and_mutate(my_gene.exons, subsequence, 19)
    return Seq(myh7_snp), Seq(mutated_results_exons[0][1])



if __name__ == "__main__":
    #combining the two latest datasets to reduce unecessery copies
    df1 = pd.read_csv('data_updated_inhibition.csv')
    df2 = pd.read_csv('data_with_flags.csv')
    extra_columns = df2.columns.difference(df1.columns)
    df_check = df1.merge(df2[['index'] + extra_columns.tolist()], on='index', how='left')

    #running find_aso function again to correct writing of missing aso
    cols_to_override = ['Transcript', 'Location_in_sequence', 'Location_div_by_length']
    df_sub = df_check.copy()
    df_sub = df_sub.drop(columns=cols_to_override, errors='ignore')

    dict_reg_id, _, hbv_data, dict_reg_seq, _ = get_relevant_data()
    df_sub_unique_aso = get_sub_df(df_sub, list(dict_reg_seq.keys()))
    df_unique_update, relevant_transcripts_dict = find_aso(df_sub_unique_aso, dict_reg_id, dict_reg_seq)
    df_sub_hbv = get_sub_df(df_sub, ['HBV'])
    df_unique_hbv = aso_for_HBV(df_sub_hbv, hbv_data)

    df_modified_temp = pd.concat([
        df_unique_update[['Sequence', 'Transcript', 'Location_in_sequence', 'Location_div_by_length']],
        df_unique_hbv[['Sequence', 'Transcript', 'Location_in_sequence', 'Location_div_by_length']]
    ], ignore_index=True)

    df_modified_temp['has_location'] = df_modified_temp['Location_in_sequence'].notna()
    df_modified_temp_sorted = df_modified_temp.sort_values(
        by=['Sequence', 'has_location'],
        ascending=[True, False]
    )
    df_modified = df_modified_temp_sorted.drop_duplicates(subset='Sequence', keep='first')
    df_modified = df_modified.drop(columns=['has_location'])
    print(f"Is 'Sequence' unique in df_modified before merge? {df_modified['Sequence'].is_unique}")

    df_with_info = df_sub.copy()
    df_with_info = df_with_info.merge(
        df_modified,
        on='Sequence',
        how='left'
    )

    check_function(df_sub, df_with_info)

    # exporting transcripts to csv
    df_transcripts = pd.DataFrame(
        relevant_transcripts_dict.items(),
        columns=['Transcript', 'Sequence']
    )
    df_transcripts.to_csv('Transcripts_for_data.csv',  index=False)

    #fixing location for known SNP in MYH7 gene
    gene_list = ['MYH7']
    df_seq = pd.read_csv('Genes with Sequences.csv')
    dict_id = {'MYH7' : df_seq[df_seq['gene'] == 'MYH7']['gene_id'].values[0]}
    dict_seq = {'MYH7' : df_seq[df_seq['gene'] == 'MYH7_SNP_715C']['gene_sequence'].values[0]}
    df_sub_unique_aso = get_sub_df(df_with_info, gene_list)
    seq_list = ['TCAGCTTGGCGATGATCT', 'CAGCTTGGCGATGATCT']
    df_sub_unique_aso = df_sub_unique_aso.loc[df_sub_unique_aso['Sequence'].isin(seq_list)]
    df_unique_update, _ = find_aso(df_sub_unique_aso, dict_id, dict_seq)
    df_unique_update.loc[df_unique_update["Sequence"] == seq_list, 'Transcript'] = 'MYH7_SNP_715C'
    df_copy = df_with_info.copy()
    for seq in seq_list:
        df_copy.loc[df_copy['Sequence'] == seq, 'Transcript'] = df_unique_update[df_unique_update['Sequence'] == seq]['Transcript'].values[0]
        df_copy.loc[df_copy['Sequence'] == seq, 'Location_in_sequence'] = df_unique_update[df_unique_update['Sequence'] == seq]['Location_in_sequence'].values[0]
        df_copy.loc[df_copy['Sequence'] == seq, 'Location_div_by_length'] = df_unique_update[df_unique_update['Sequence'] == seq]['Location_div_by_length'].values[0]


    # Compare relevant columns
    mask_not_in_seq = ~df_with_info['Sequence'].isin(seq_list)
    unchanged = df_with_info.loc[mask_not_in_seq, ['Transcript', 'Location_in_sequence', 'Location_div_by_length']].equals(
        df_copy.loc[mask_not_in_seq, ['Transcript', 'Location_in_sequence', 'Location_div_by_length']]
    )
    if unchanged:
        print("✅ All rows NOT in seq_list remained unchanged.")
    else:
        print("❌ Some rows NOT in seq_list were changed.")

    checking = df_copy[df_copy['Sequence'].isin(seq_list)].sample(10, random_state=42)[['Sequence', 'Transcript', 'Location_in_sequence', 'Location_div_by_length']]


    df_copy.to_csv('data_asoptimizer_05_08.csv', index=False)




