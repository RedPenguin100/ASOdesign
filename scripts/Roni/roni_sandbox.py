import pandas as pd
import numpy as np
from pandas.core.interchange.dataframe_protocol import DataFrame
import pickle
from off_target_functions import parse_fasta, parse_gtf
from Bio import SeqIO



path_1 = "/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/"
path_2 = "/home/oni/iGEM/data/"
# print('test')
# # parse_gtf(path_1 +'gencode.v48.chr_patch_hapl_scaff.annotation.gtf',
# #           path_2+'_gtf_annotations.pkl')
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

df = pd.read_csv(path_1+'ACH-000232_transcriptome_premRNA.merged.csv')
print(df.shape)