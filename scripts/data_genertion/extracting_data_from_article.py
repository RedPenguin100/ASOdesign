import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def load_csv(file_path):
    """
    Reads a CSV file and returns a pandas DataFrame.
    :param file_path: Path to the CSV file
    :return: DataFrame containing CSV data
    """
    try:
        df = pd.read_csv(file_path)
        print("CSV file loaded successfully!")
        return df
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return None

def filter_by_mod(df, mod):
    # Observing one ASO and the same modification, in different locations
    filtered_df = df[df["Modification"] == mod]
    sorted_filtered_df = filtered_df.sort_values(by="Inhibition(%)", ascending=False)
    return sorted_filtered_df

def create_summary(df):
    # Analysing the data to understand how many unique values there are in genes, sequences and chemical modifications
    unique_target_genes = df['Target_gene'].unique()
    unique_sequence = df['Sequence'].unique()
    unique_mod = df['Modification'].unique()
    unique_linkage = df['Linkage'].unique()
    unique_cell_line = df['Cell_line'].unique()
    cell_line_counts = df['Cell_line'].value_counts().reindex(unique_cell_line, fill_value=0)
    # creating new df as summary of the columns and its meaning
    meaning_and_examples = ['Serial number based on a company named ISIS', ', '.join(unique_target_genes.tolist()),
                            'Cell lines', 'Density', 'Transfection', 'Volume', 'Treatment period', 'Primer probe set',
                            'ASO sequences (15470 unique sequences)', ', '.join(unique_mod.tolist()),
                            'Modification Locations (1561 unique locations)', 'Chemical Pattern',
                            ', '.join(unique_linkage.tolist()),
                            'Location of linkages (47 unique locations)',
                            'String that represents the 3D structure of sequence',
                            'Inhibition level measured by real time qRT-PCR', 'Sequence length']
    df_summary = pd.DataFrame({
        'Index': range(len(df.dtypes)),
        'Column Name': df.columns,
        'Data Type': df.dtypes.values,
        'Explanation': meaning_and_examples
    })
    seq_counts = df['Sequence'].value_counts().reindex(unique_sequence, fill_value=0)  # counting repetitive sequences
    df_aso_seq = pd.DataFrame({'ASO Sequences': unique_sequence,
                               'Number of repetitions': seq_counts})  # creating a df to summaries repetitions
    aso_count_to_gene = df.groupby('Target_gene')['Sequence'].nunique().reindex(unique_target_genes,
                                                                                fill_value=0)  # counting how many unique aso a gene has
    df_genes = pd.DataFrame({'Genes': unique_target_genes,
                             'Number of ASO per gene': aso_count_to_gene})  # creating a df of unique aso count to gene
    df_cell_line = pd.DataFrame({'Cell Lines': unique_cell_line, 'Count': cell_line_counts})  # creating a df of unique cell lines
    # exporting excal file of summary of columns and analysis data
    output_file = 'Summary of DataBase.xlsx'
    with pd.ExcelWriter(output_file, engine="xlsxwriter") as writer:
        df_summary.to_excel(writer, sheet_name="Total Summary", index=False)
        df_genes.to_excel(writer, sheet_name="Genes List", index=False)
        df_aso_seq.to_excel(writer, sheet_name="ASOs Sequences", index=False)
        df_cell_line.to_excel(writer, sheet_name="Cell Lines", index=False)


def find_if_seq_in_gene(seq, gene):
    # receives a sequence and a gene
    # if sequence in gene, returns the index, else returns -1
    seq_up = seq.upper()
    gene_up = gene.upper()
    return gene_up.find(seq_up)

def read_fasta_biopython(file_path):
    # reads fasta file and returns dictionary of sequence and its name (record.id)
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def filter_cell_line_human(df, cell_line_info):
    # receives a df of all data and a df of the cell lines and the organism their belong to
    # returns only the rows with human cell lines
    cell_human = cell_line_info.loc[cell_line_info['Cell line organism'] == 'human', 'Cell_line'].tolist()
    return df[df['Cell_line'].isin(cell_human)]


def which_gene_pkk(df):
    df_gene_test = df[df['Target_gene'] == 'PKK']
    pkk_aso = df_gene_test['Sequence'].tolist()

    ripk4 = read_fasta_biopython('sequence_ripk4.fasta')
    gene_ripk4 = list(ripk4.values())[0]
    klkb1 = read_fasta_biopython('sequence_klkb1.fasta')
    gene_klkb1 = list(klkb1.values())[0]
    index_ripk4 = []
    index_klkb1 = []
    for aso in pkk_aso:
        this_aso = Seq(aso)
        index_ripk4.append(find_if_seq_in_gene(str(this_aso.reverse_complement()), gene_ripk4))
        index_klkb1.append(find_if_seq_in_gene(str(this_aso.reverse_complement()), gene_klkb1))
    result_ripk4 = any(num != -1 for num in index_ripk4)
    result_klkb1 = any(num != -1 for num in index_klkb1)
    return {'RIPK4': result_ripk4, 'KLKB1': result_klkb1 }


def which_gene_tau(df):
    df_gene_test = df[df['Target_gene'] == 'Tau']
    tau_aso = df_gene_test['Sequence'].tolist()

    mapt = read_fasta_biopython('sequence_MAPT.fasta')
    gene_mapt = list(mapt.values())[0]
    index_mapt = []
    for aso in tau_aso:
        this_aso = Seq(aso)
        index_mapt.append(find_if_seq_in_gene(str(this_aso.reverse_complement()), gene_mapt))
    result_mapt = any(num != -1 for num in index_mapt)
    return {'MAPT': result_mapt, 'index':index_mapt}


def gene_snca(df):
    df_gene_test = df[df['Target_gene'] == 'SNCA']
    snca_aso = df_gene_test['Sequence'].tolist()

    snca = read_fasta_biopython('sequence_snca.fasta')
    gene_snca = list(snca.values())[0]
    index_snca = []
    for aso in snca_aso:
        this_aso = Seq(aso)
        index_snca.append(find_if_seq_in_gene(str(this_aso.reverse_complement()), gene_snca))
    result_snca = any(num != -1 for num in index_snca)
    return result_snca



if __name__ == "__main__":
    file_path = "experiments_with_smiles.csv"
    df = load_csv(file_path)
    if df is not None:
        print(df.head())  # Display the first few rows
        print(df.dtypes)

        # When loading the DF at first, an error occurred stating some columns had mixed data types in them
        # Identifying the columns with mixed data types and converting for one type
        mixed_cols = [col for col in df.columns if df[col].map(type).nunique() > 1]
        print("Columns with mixed data types:", mixed_cols)
        for col in mixed_cols:
            print(f"Column: {col}")
            print(df[col].unique()[:10])  # Show first 10 unique values
        df[mixed_cols[0]] = pd.to_numeric(df[mixed_cols[0]], errors='coerce')
        df[mixed_cols[1]] = df[mixed_cols[1]].astype(str)
        df[mixed_cols[2]] = df[mixed_cols[2]].astype(str)
        mixed_cols_new = [col for col in df.columns if df[col].map(type).nunique() > 1]
        print("Columns with mixed data types:", mixed_cols_new)
        print(df.dtypes)
        create_summary(df)

        df = df.reset_index()
        #testing search for gene functions with known gene
        test_func_snca = gene_snca(df)

        # adding a column of the canonical name of the genes
        # if original name is inconclusive, a check will run and determine based on ASOs sequences in df
        where_is_pkk = which_gene_pkk(df)
        is_tau_in_mapt = which_gene_tau(df)
        df_genes_info = load_csv('genes_info.csv')
        df_new = df.merge(df_genes_info[['Target_gene', 'Canonical Gene Name']], on='Target_gene', how='left') # note that 2 genes were determined by their ASOs sequences, and 1 with the help of the internet

        # to add a column of cell line organism to original df
        # some of the cell lines' organism are inconclusive
        df_cell_line_info = load_csv('cell lines info.csv')
        merged_df = df_new.merge(df_cell_line_info, on='Cell_line', how='left')

        # to filter df only for human cell lines, run the commented lines
        # some of the cell lines' organism are inconclusive
        df_filtered = filter_cell_line_human(df_new, df_cell_line_info)

        merged_df.to_csv("data_from_article_fixed.csv", index=False)
