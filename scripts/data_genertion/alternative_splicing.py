import gget
from scripts.data_genertion.extracting_data_from_article import find_if_seq_in_gene, load_csv, read_fasta_biopython
from Bio.Seq import Seq
import gffutils
import os
import pandas as pd



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
            in_this_transcript = find_if_seq_in_gene(str(reverse_aso_seq), gene_transcripts[i])
            if in_this_transcript != -1:
                return in_this_transcript, gene_transcripts[i-1], gene_transcripts[i]
    return -1, -1, -1

def make_transcripts_dict(locus_to_data):
    gene_transcripts = {}
    for gene_name in locus_to_data:
        #gene_transcripts[gene_name] = get_transcripts_of_gene(locus_to_data[gene_name].gene_id)
        gene_transcripts[gene_name] = get_transcripts_of_gene(locus_to_data[gene_name])
    return gene_transcripts


def find_aso(df, locus_to_data, df_with_seq):
    # receiving dataframe with sequences of aso and their canonical gene name, dictionary with gene info and df of genes with their sequence
    # returning the dataframe with added columns about the first location of an aso
    # only for human sequences
    gene_transcripts_dict = make_transcripts_dict(locus_to_data)
    aso_unique_sequence = df['Sequence'].unique()
    for aso in aso_unique_sequence:
        reverse_aso_seq = str(Seq(aso).reverse_complement())
        this_rows = df[df['Sequence'] == aso]
        this_aso_gene = this_rows["Canonical Gene Name"].iloc[0]
        if this_aso_gene in locus_to_data:
            gene_seq = str(df_with_seq[df_with_seq["gene"] == this_aso_gene]['gene_sequence'])
            in_pre_mrna = find_if_seq_in_gene(reverse_aso_seq, gene_seq)
            if in_pre_mrna != -1:
               df.loc[df["Sequence"] == aso, 'Where_ASO_is_found'] = 'pre_mRNA'
               df.loc[df["Sequence"] == aso, 'Seqeuence_of_location'] = gene_seq
               df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = int(in_pre_mrna)
            else:
                in_transcript, transcript_id, transcript_seq = find_in_transcripts(reverse_aso_seq, this_aso_gene,
                                                                                   gene_transcripts_dict)
                if in_transcript != -1:
                    df.loc[df["Sequence"] == aso, 'Where_ASO_is_found'] = transcript_id
                    df.loc[df["Sequence"] == aso, 'Seqeuence_of_location'] = transcript_seq
                    df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = int(in_transcript)
        else:
            df.loc[df["Sequence"] == aso, 'Where_ASO_is_found'] = None
            df.loc[df["Sequence"] == aso, 'Seqeuence_of_location'] = None
            df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = None
    return df

def find_aso_for_snord(df, locus_to_data, df_with_seq):
    gene_transcripts_dict = make_transcripts_dict(locus_to_data)
    aso_unique_sequence = df['Sequence'].unique()
    for aso in aso_unique_sequence:
        reverse_aso_seq = str(Seq(aso).reverse_complement())
        this_rows = df[df['Sequence'] == aso]
        this_aso_gene = this_rows["Canonical Gene Name"].iloc[0]
        if this_aso_gene in locus_to_data:
            gene_seq = str(df_with_seq[df_with_seq["gene"] == this_aso_gene]['gene_sequence'])
            in_pre_mrna = find_if_seq_in_gene(reverse_aso_seq, gene_seq)
            if in_pre_mrna != -1:
                df.loc[df["Sequence"] == aso, 'Where_ASO_is_found'] = 'pre_mRNA'
                df.loc[df["Sequence"] == aso, 'Seqeuence_of_location'] = gene_seq
                df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = int(in_pre_mrna)
            else:
                in_transcript, transcript_id, transcript_seq = find_in_transcripts(reverse_aso_seq, this_aso_gene,
                                                                                   gene_transcripts_dict)
                if in_transcript != -1:
                    df.loc[df["Sequence"] == aso, 'Where_ASO_is_found'] = transcript_id
                    df.loc[df["Sequence"] == aso, 'Seqeuence_of_location'] = transcript_seq
                    df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = int(in_transcript)
        else:
            df.loc[df["Sequence"] == aso, 'Where_ASO_is_found'] = None
            df.loc[df["Sequence"] == aso, 'Seqeuence_of_location'] = None
            df.loc[df["Sequence"] == aso, 'Location_in_sequence'] = None
    return df

def check_get_id_from_annotations(file_path, db_path, gene_list):
    if not os.path.exists(db_path):
        print("Creating DB...")
        gffutils.create_db(
            file_path,
            dbfn=db_path,
            force=True,
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
            merge_strategy='merge',
            verbose=True,
        )
        print("DB created")
    else:
        print("DB already exists. Using cached version.")

    db = gffutils.FeatureDB(db_path)
    gene_name_to_id = {}

    for gene in db.features_of_type('gene'):
        # Safely extract gene_id and gene_name
        gene_id = gene.attributes.get('gene_id', [None])[0]
        gene_name = gene.attributes.get('gene_name', [None])[0]

        # Only add if both are present
        if gene_id and gene_name and gene_name in gene_list:
            gene_name_to_id[gene_name] = gene_id
    return gene_name_to_id

def get_seq_for_genes(dict_gene):
    df = pd.DataFrame(list(dict_gene.items()), columns=["gene", "gene_id"])
    df['gene_sequence'] = None
    for gene in dict_gene:
        seq_get = gget.seq(dict_gene[gene])
        seq = seq_get[1]
        df.loc[df["gene"] == gene, 'gene_sequence'] = seq
    return df


if __name__ == "__main__":
    df = load_csv('data_from_article_fixed.csv')
    gene_list = df['Canonical Gene Name'].unique()
    locus_to_data = check_get_id_from_annotations('gencode.v34.basic.annotation.gtf', 'human_annotations.db', gene_list)
    df_with_seq = get_seq_for_genes(locus_to_data)
    snord115_info = gget.search('SNORD115', species='homo_sapiens')
    snord_dict = snord115_info.set_index('gene_name')['ensembl_id'].to_dict()
    df_snord_with_seq = get_seq_for_genes(snord_dict)

    # df_sub_unique_aso = get_sub_df(df, df_with_seq['gene'].tolist())
    # df_unique_update = find_aso(df_sub_unique_aso, locus_to_data, df_with_seq)
    # df_sub_snord = get_sub_df(df, 'HBII-52')
    #

    HBV = read_fasta_biopython('GCF_HBV.fna')
    seq_HBV = list(HBV.values())[0]
    df_with_seq.loc[len(df_with_seq)] = ['HBV', None, seq_HBV]
    df_genes_seq = pd.concat([df_with_seq, df_snord_with_seq], ignore_index=True)
    df_genes_seq.to_csv('Genes with Sequences.csv', index=False)
