import pandas as pd
import numpy as np
from pandas.core.interchange.dataframe_protocol import DataFrame
import pickle

#from ASOdesign.scripts.Roni.mutate_cell_line_transcriptome import celline_list
from off_target_functions import parse_fasta, parse_gtf
from Bio import SeqIO
from itertools import islice

#from scripts.Roni.gtf.compare_gtf import path_1

path_1 = "/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/"
path_2 = "/home/oni/ASOdesign/scripts/Roni/"
path_3 = "/home/oni/ASOdesign/scripts/Roni/gtf/"

print('run01')
# with open(path_1+'_gtf_annotations.pkl', 'rb') as f:
#     gtf_dict = pickle.load(f)
# print(gtf_dict['ENST00000340650.6']) # 71 out of range
# print(gtf_dict['ENST00000367976.4']) # 314 out of range
# print(gtf_dict['ENST00000164024.5']) # 7648 out of range
# print(gtf_dict['ENST00000373345.9']) # 43438 out of range
# print(gtf_dict['ENST00000265293.9']) # 105646 same shift worked on the first mutation

# df = pd.read_csv(path_1+'ACH-001328.mutated_transcriptome_premRNA.merged.csv')
# print(df[df['Transcript_ID'] == 'ENST00000373345.9']['Original Transcript Sequence'].values[0])
# print(df[df['Transcript_ID'] == 'ENST00000373345.9']['Mutated Transcript Sequence'].values[0])
#
# print(df['Mutated Transcript Sequence'].count())

#======================================================================================================================
#======================================================================================================================
#================================== PARSE GTF ================================================================

print('test')
parse_gtf(path_3 +'gencode.v34.chr_patch_hapl_scaff.annotation.gtf',
          path_2+'_gtf_annotations_v34.pkl')

#======================================================================================================================
#======================================================================================================================
#
# with open(path_1+'_whole_genomic_sequence.pkl', 'rb') as f:
#     genome_dict = pickle.load(f)
# print('Read GTF.')
#
# excluded_chroms = {f'CM{str(i).zfill(6)}.2' for i in range(663, 687)}
# small_dict = {k: v for k, v in genome_dict.items() if k not in excluded_chroms}
#
# with open(path_1+'chr_M_rest.pkl', "wb") as f:
#     pickle.dump(small_dict, f)
#
#
# #parse_fasta(path_2+'GCA_000001405.29_GRCh38.p14_genomic.fna', path_1+'_whole_genomic_sequence.pkl')
# #parse_fasta(path_2+'Homo_sapiens.GRCh38.dna_sm.toplevel.fa', path_1+'_full_genomic_sequence.pkl')
# print('worked')#test
# print()


#
# lennn = gtf_dict['ENST00000396024.7']['end'] - gtf_dict['ENST00000396024.7']['start']
# print(f'real gene len is {lennn+1}')
# df = pd.read_csv(path_1+'ACH-000232_transcriptome_premRNA.merged.csv')
# row = df[df['Transcript_ID'] == 'ENST00000396024.7']
#
# # Extract the sequence
# sequence = row['Original Transcript Sequence'].values[0]
#
# # Optionally print its length
# print(f"Sequence length: {len(sequence)}")
#
# print()


# with open(path_2+'off_target_info_mRNA.pkl', 'rb') as f:
#     genome_dict = pickle.load(f)
#
# genome_dict[6].to_csv('test.csv')

#======================================================================================================================
#======================================================================================================================
#================================== CREATE SCATTER PLOT ================================================================

# from asodesigner.consts import DATA_PATH_NEW
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # Load both CSVs
# df1 = pd.read_csv(DATA_PATH_NEW/"off_target.GEO.top100.cutoff1000.premRNA.csv", index_col=0)  # set index to the column with gene IDs
# df2 = pd.read_csv(DATA_PATH_NEW/"data_asoptimizer_updated.csv", index_col=0)
#
# # Pick the columns you want to compare
# col1 = "off_target_score_MS"
# col2 = "Inhibition(%)"
#
# # Merge on the index
# merged = df1[[col1]].join(df2[[col2]], how="inner")
#
# # Scatter plot
# plt.figure(figsize=(6,6))
# plt.scatter(merged[col1], merged[col2], alpha=0.5)
# plt.xlabel(col1)
# plt.ylabel(col2)
# plt.title(f"Scatter plot: {col1} vs {col2}")
# plt.grid(True)
# plt.show()

#======================================================================================================================
#======================================================================================================================

# with open('/home/oni/ASOdesign/scripts/data_genertion/off_target/off_target_info.premRNA.top100.cutoff1200.pkl', 'rb') as f:
#      dicto = pickle.load(f)
#
# pd.set_option('display.max_columns', None)
#
# # optionally show all rows (careful if DataFrame is huge)
# pd.set_option('display.max_rows', None)
#
# # widen the display so it doesn't wrap badly
# pd.set_option('display.width', 200)
#
# #print(dicto)
#
# print(len(dicto.keys()))

#======================================================================================================================
#======================================================================================================================
#================================== COMPARE SEQUENCES  ================================================================

celline_list = ["ACH-000232", "ACH-000463", "ACH-000739", "ACH-001086", "ACH-001188", "ACH-001328"]

for cell in celline_list:
    count = 0
    path34 = path_2 + cell + "_transcriptome_premRNA.merged.csv"
    v34 = pd.read_csv(path34)
    print("read v34")
    path48 = path_1 + cell + ".mutated_transcriptome_premRNA.merged.csv"
    v48 = pd.read_csv(path48)
    print("read v48")

    # make Transcript_ID the index for fast lookup
    v34 = v34.set_index("Transcript_ID")

    for idx, row in v48.iterrows():
        t_id = row["Transcript_ID"]
        seq48 = row["Original Transcript Sequence"]

        # get seq34 if transcript exists
        if t_id not in v34.index:
            continue
        seq34 = v34.at[t_id, "Original Transcript Sequence"]

        # skip if either sequence is NaN
        if pd.isna(seq34) or pd.isna(seq48):
            continue

        # check mismatch
        if seq34 != seq48:
            count += 1
            print(f"Transcript: {t_id}\n"
                  f"seq34:\n s:{seq34[:100]}\n e:{seq34[-100:]}\n len:{len(seq34)}\n\n"
                  f"seq48:\n s:{seq48[:100]}\n e:{seq48[-100:]}\n len:{len(seq48)}\n\n")

    print(f"cell line: {cell} has {count} mismatches")



