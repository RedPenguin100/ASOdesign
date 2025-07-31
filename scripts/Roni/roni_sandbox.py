import pandas as pd
import numpy as np
from pandas.core.interchange.dataframe_protocol import DataFrame
import pickle
from off_target_functions import parse_fasta, parse_gtf
from Bio import SeqIO



path_1 = "/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/"
path_2 = "/home/oni/iGEM/data/"

# #parse_gtf(path +'gencode.v48.chr_patch_hapl_scaff.annotation.gtf',
# #           '../data_genertion/cell_line_expression/_NEW_gtf_annotations.pkl')
#
# with open(path+'_NEW_gtf_annotations.pkl', 'rb') as f:
#     gtf_dict = pickle.load(f)
# print('Read GTF.')
#
#
# print(gtf_dict['ENST00000474973.5'])
print('test')
#parse_fasta(path_2+'GCA_000001405.29_GRCh38.p14_genomic.fna', path_1+'_whole_genomic_sequence.pkl')
parse_fasta(path_2+'Homo_sapiens.GRCh38.dna_sm.toplevel.fa', path_1+'_full_genomic_sequence.pkl')
print('worked')