import numpy as np
from external.risearch.RIsearch1.numpy_to_csv import dsm_variable_to_csv,numpy_to_csv

# got error when did option 2 np_array.transpose(3, 2, 1, 0)
def rna_dna_to_dna_rna(np_array):
    return np_array.transpose(2, 3, 0, 1)
def rna_dna_to_ps_dna_rna(np_array,eps_array):
    return  rna_dna_to_dna_rna(np_array) + eps_array

def make_dsm_dna_rna():
    dsm_variable_to_csv(dsm_name = "dsm_su95_rev_woGU_pos",func=rna_dna_to_dna_rna, kwargs=None, func_extend=rna_dna_to_dna_rna, kwargs_extend=None)
def make_dsm_ps_dna_rna(eps_array =  np.full((6, 6, 6, 6), 1.3) ):
        dsm_variable_to_csv(dsm_name = "dsm_su95_rev_woGU_pos",func=rna_dna_to_ps_dna_rna, kwargs={'eps_array' : eps_array} , func_extend=rna_dna_to_dna_rna, kwargs_extend=None)

