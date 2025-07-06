import gget
from scripts.data_genertion.extracting_data_from_article import find_if_seq_in_gene, load_csv, read_fasta_biopython
from Bio.Seq import Seq
import gffutils
import os
import pandas as pd
from asodesigner.read_human_genome import get_locus_to_data_dict



def get_transcripts_of_gene(gene_ID):
    # using gene Ensembl ID to get gene transcripts as a list
    return gget.seq(gene_ID, isoforms=True)

def get_sub_df(df, gene_list):
    # receiving df, and returning sub-df with only unique ASO and their canonical gene names
    filtered_df = df[df["Canonical Gene Name"].isin(gene_list)]
    unique_df = filtered_df.drop_duplicates(subset="Sequence", keep="first")
    return unique_df[["Sequence", "Canonical Gene Name"]].reset_index(drop=True)

def find_in_transcripts(reverse_aso_seq, gene_name, transcripts_dict):
    if gene_name in transcripts_dict:
        gene_transcripts = transcripts_dict[gene_name]
        for i in range(1,len(gene_transcripts),2):
            transcript_seq = gene_transcripts[i]
            in_this_transcript = find_if_seq_in_gene(str(reverse_aso_seq), transcript_seq)
            if in_this_transcript != -1:
                gene_id = gene_transcripts[i-1]
                gene_id = gene_id[1:].split()[0]
                return in_this_transcript, gene_id, transcript_seq
    return -1, -1, -1

def make_transcripts_dict(locus_to_data):
    gene_transcripts = {}
    for gene_name in locus_to_data:
        gene_transcripts[gene_name] = get_transcripts_of_gene(locus_to_data[gene_name])
    return gene_transcripts


def find_aso(df, locus_to_data, dict_with_seq):
    # receiving dataframe with sequences of aso and their canonical gene name, dictionary with gene info and df of genes with their sequence
    # returning the dataframe with added columns about the first location of an aso
    # only for human sequences
    gene_transcripts_dict = make_transcripts_dict(locus_to_data)
    aso_unique_sequence = df['Sequence'].unique()
    relevant_transcripts_dict = {}
    for aso in aso_unique_sequence:
        reverse_aso_seq = str(Seq(aso).reverse_complement())
        this_rows = df[df['Sequence'] == aso]
        this_aso_gene = this_rows["Canonical Gene Name"].iloc[0]
        found = False
        if this_aso_gene in locus_to_data:
            gene_seq = str(dict_with_seq[this_aso_gene])
            in_pre_mrna = find_if_seq_in_gene(reverse_aso_seq, gene_seq)
            if in_pre_mrna != -1:
               df.loc[df["Sequence"] == aso, 'Transcript'] = None
               df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = int(in_pre_mrna)
               df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = int(in_pre_mrna)/len(gene_seq)
               found = True
            else:
                in_transcript, transcript_id, transcript_seq = find_in_transcripts(reverse_aso_seq, this_aso_gene,
                                                                                   gene_transcripts_dict)
                if in_transcript != -1:
                    df.loc[df["Sequence"] == aso, 'Transcript'] = transcript_id
                    df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = int(in_transcript)
                    df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = int(in_transcript)/len(transcript_seq)
                    relevant_transcripts_dict[transcript_id] = transcript_seq
                    found = True
                    # all transcripts that exist, create a different csv with transcript id and its sequence, only relevant ones
        if not found:
            df.loc[df["Sequence"] == aso, 'Transcript'] = None
            df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = None
            df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = None
    return df, relevant_transcripts_dict


def find_aso_for_snord(df, locus_to_data, dict_with_seq):
    relevant_transcripts_dict = {}
    aso_unique_sequences = df['Sequence'].unique()

    for aso in aso_unique_sequences:
        reverse_aso_seq = str(Seq(aso).reverse_complement())
        found = False
        for key, value in dict_with_seq.items():
            gene_seq = str(value)
            this_id = locus_to_data.get(key)
            # First try to find in gene sequence
            pos_in_gene = find_if_seq_in_gene(reverse_aso_seq, gene_seq)
            if pos_in_gene != -1:
                df.loc[df["Sequence"] == aso, 'Transcript'] = this_id
                df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = pos_in_gene
                df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = pos_in_gene / len(gene_seq)
                found = True
                break
            # Try to find in transcripts
            gene_transcripts_dict = make_transcripts_dict({key: this_id})
            pos_in_transcript, transcript_id, transcript_seq = find_in_transcripts(reverse_aso_seq, key, gene_transcripts_dict)
            if pos_in_transcript != -1:
                df.loc[df["Sequence"] == aso, 'Transcript'] = transcript_id
                df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = pos_in_transcript
                df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = pos_in_transcript / len(transcript_seq)
                relevant_transcripts_dict[transcript_id] = transcript_seq
                found = True
                break
        # If not found in any gene or transcript
        if not found:
            df.loc[df["Sequence"] == aso, 'Transcript'] = None
            df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = None
            df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = None
    return df, relevant_transcripts_dict

def aso_for_HBV(df, hbv_data):
    aso_unique_sequences = df['Sequence'].unique()
    for aso in aso_unique_sequences:
        reverse_aso_seq = str(Seq(aso).reverse_complement())
        hbv_seq = hbv_data.iloc[0]['gene_sequence']
        where_aso = find_if_seq_in_gene(reverse_aso_seq, hbv_seq)
        if where_aso != -1:
            df.loc[df["Sequence"] == aso, 'Transcript'] = None
            df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = int(where_aso)
            df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = int(where_aso) / len(hbv_seq)
        else:
            df.loc[df["Sequence"] == aso, 'Transcript'] = None
            df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = None
            df.loc[df["Sequence"] == aso, 'Location_div_by_length'] = None
    return df



def get_id_from_annotations(file_path, gene_list):
    gene_name_to_id = {}

    with open(file_path, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue  # skip header

            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue  # skip malformed lines

            feature_type = columns[2]
            if feature_type != 'gene':
                continue

            # Parse the attributes field
            attributes_field = columns[8]
            attributes = {}
            for attr in attributes_field.split(';'):
                if attr.strip():
                    key_value = attr.strip().split(' ', 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        attributes[key] = value.strip('"')

            gene_name = attributes.get('gene_name')
            gene_id = attributes.get('gene_id')

            if gene_name and gene_id and gene_name in gene_list:
                gene_name_to_id[gene_name] = gene_id

    return gene_name_to_id

def get_seq_for_genes(dict_with_id, locus_to_data):
    df = pd.DataFrame(list(dict_with_id.items()), columns=["gene", "gene_id"])
    df['gene_sequence'] = None
    for gene in dict_with_id:
        seq = str(locus_to_data[gene].full_mrna)
        df.loc[df["gene"] == gene, 'gene_sequence'] = seq
    return df

def get_seq_for_snord(dict_gene):
    df = pd.DataFrame(list(dict_gene.items()), columns=["gene", "gene_id"])
    df['gene_sequence'] = None
    for gene in dict_gene:
        seq_get = gget.seq(dict_gene[gene])
        seq = seq_get[1]
        df.loc[df["gene"] == gene, 'gene_sequence'] = seq
    return df


def get_initial_sequences(df, gene_list):
    # run this only one time to load gene sequences to csv, the file is already exist after
    # for this function to work, need to upload to this directory 'gencode.v34.basic.annotation.gtf'
    id_to_data = get_id_from_annotations('gencode.v34.basic.annotation.gtf', gene_list)
    locus_to_data = get_locus_to_data_dict(gene_subset=gene_list)
    df_with_seq = get_seq_for_genes(id_to_data, locus_to_data)
    snord115_info = gget.search('SNORD115', species='homo_sapiens')
    snord_dict = snord115_info.set_index('gene_name')['ensembl_id'].to_dict()
    df_snord_with_seq = get_seq_for_snord(snord_dict)
    HBV = read_fasta_biopython('GCF_HBV.fna')
    seq_HBV = list(HBV.values())[0]
    df_with_seq.loc[len(df_with_seq)] = ['HBV', None, seq_HBV]
    df_genes_seq = pd.concat([df_with_seq, df_snord_with_seq], ignore_index=True)
    df_genes_seq.to_csv('Genes with Sequences.csv', index=False)


def get_relevant_data():
    df_genes_seq = load_csv('Genes with Sequences.csv')
    df_reg_genes = df_genes_seq[:18]
    hbv_data = df_genes_seq[18:19]
    df_snords = df_genes_seq[19:]
    dict_reg_id = dict(zip(df_reg_genes['gene'], df_reg_genes['gene_id']))
    dict_snord_id = dict(zip(df_snords['gene'], df_snords['gene_id']))
    dict_reg_seq = dict(zip(df_reg_genes['gene'], df_reg_genes['gene_sequence']))
    dict_snord_seq = dict(zip(df_snords['gene'], df_snords['gene_sequence']))
    return dict_reg_id, dict_snord_id, hbv_data, dict_reg_seq, dict_snord_seq




if __name__ == "__main__":
    df = load_csv('data_from_article_fixed.csv')
    gene_list = df['Canonical Gene Name'].unique()
    #get_initial_sequences(df, gene_list) ## for this function to work, need to upload to this directory 'gencode.v34.basic.annotation.gtf'
    df_genes_seq_n = load_csv('Genes with Sequences.csv')
    dict_reg_id, dict_snord_id, hbv_data, dict_reg_seq, dict_snord_seq = get_relevant_data()
    df_sub_unique_aso = get_sub_df(df, list(dict_reg_seq.keys()))
    df_unique_update, relevant_transcripts_dict = find_aso(df_sub_unique_aso, dict_reg_id, dict_reg_seq)
    df_sub_snord = get_sub_df(df, ['HBII-52'])
    df_unique_snord, transcripts_snord = find_aso_for_snord(df_sub_snord, dict_snord_id, dict_snord_seq)
    df_sub_hbv = get_sub_df(df, ['HBV'])
    df_unique_hbv = aso_for_HBV(df_sub_hbv, hbv_data)

    df_modified_temp = pd.concat([
        df_unique_update[['Sequence', 'Transcript', 'Location_in_sequence', 'Location_div_by_length']],
        df_unique_snord[['Sequence', 'Transcript', 'Location_in_sequence', 'Location_div_by_length']],
        df_unique_hbv[['Sequence', 'Transcript', 'Location_in_sequence', 'Location_div_by_length']]
    ], ignore_index=True)

    # --- NEW LOGIC TO PRIORITIZE NON-NAN LOCATION ---
    # 1. Create a helper column that is False for NaN and True for non-NaN
    df_modified_temp['has_location'] = df_modified_temp['Location_in_sequence'].notna()

    # 2. Sort: True (non-NaN) comes before False (NaN).
    # Sorting by 'Sequence' secondarily ensures stable sorting for identical sequences.
    df_modified_temp_sorted = df_modified_temp.sort_values(
        by=['Sequence', 'has_location'],
        ascending=[True, False]  # Sort 'Sequence' ascending, 'has_location' descending (True first)
    )

    # 3. Drop duplicates: keep='first' will now keep the row with non-NaN location if one exists.
    df_modified = df_modified_temp_sorted.drop_duplicates(subset='Sequence', keep='first')

    # 4. Drop the temporary 'has_location' column
    df_modified = df_modified.drop(columns=['has_location'])
    # --- END NEW LOGIC ---

    print(f"Is 'Sequence' unique in df_modified before merge? {df_modified['Sequence'].is_unique}")

    df_with_info = df.copy()
    df_with_info = df_with_info.merge(
        df_modified,
        on='Sequence',
        how='left'
    )

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


    df_with_info.to_csv('data_updated_18.5.csv', index=False)



