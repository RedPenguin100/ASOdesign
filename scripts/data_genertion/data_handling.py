import primer3
from ViennaRNA import RNA

from scripts.data_genertion.consts import CELL_LINE_ORGANISM, INHIBITION, CANONICAL_GENE
import numpy as np
import pandas as pd
from scripts.data_genertion.consts import *


def get_unique_human_genes(all_data):
    all_data_human = all_data[all_data[CELL_LINE_ORGANISM] == 'human']
    all_data_human_no_nan = all_data_human.dropna(subset=[INHIBITION]).copy()

    genes = all_data_human_no_nan[CANONICAL_GENE].copy()
    genes_u = list(set(genes))

    genes_u.remove('HBV')
    genes_u.remove('negative_control')

    return genes_u


def get_gene_to_data(genes_u):
    from asodesigner.read_human_genome import get_locus_to_data_dict
    import pickle
    from asodesigner.consts import CACHE_DIR

    cache_path = CACHE_DIR / 'gene_to_data_simple_cache.pickle'

    # TODO: hash the pickled file to avoid mis-reads
    if not cache_path.exists():
        gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=genes_u)
        with open(cache_path, 'wb') as f:
            pickle.dump(gene_to_data, f)
    else:
        with open(cache_path, 'rb') as f:
            gene_to_data = pickle.load(f)
    return gene_to_data


def get_populated_df_with_structure_features(df, genes_u, gene_to_data):
    """
    Populate "the data" df with features like exon/intron, start of the sense strand, if found.
    """
    from asodesigner.util import get_antisense
    df_copy = df.copy()
    all_data_human = df_copy[df_copy[CELL_LINE_ORGANISM] == 'human']
    all_data_human_no_nan = all_data_human.dropna(subset=[INHIBITION]).copy()
    all_data_human_gene = all_data_human_no_nan[all_data_human_no_nan[CANONICAL_GENE].isin(genes_u)].copy()
    SENSE_START = 'sense_start'
    SENSE_LENGTH = 'sense_length'
    SENSE_TYPE = 'sense_type'

    found = 0
    not_found = 0
    all_data_human_gene[SENSE_START] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_LENGTH] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_TYPE] = "NA"
    for index, row in all_data_human_gene.iterrows():
        gene_name = row[CANONICAL_GENE]
        locus_info = gene_to_data[gene_name]
        pre_mrna = locus_info.full_mrna
        antisense = row[SEQUENCE]
        sense = get_antisense(antisense)
        idx = pre_mrna.find(sense)
        all_data_human_gene.loc[index, SENSE_START] = idx
        all_data_human_gene.loc[index, SENSE_LENGTH] = len(antisense)
        if idx != -1:
            genome_corrected_index = idx + locus_info.exon_indices[0][0]
            found = False
            for exon_indices in locus_info.exon_indices:
                # print(exon[0], exon[1])
                if exon_indices[0] <= genome_corrected_index <= exon_indices[1]:
                    all_data_human_gene.loc[index, SENSE_TYPE] = 'exon'
                    found = True
                    break
        if not found:
            all_data_human_gene.loc[index, SENSE_TYPE] = 'intron'
    return all_data_human_gene


def get_populate_fold(df, genes_u, gene_to_data, fold_variants=[(40, 15)]):
    from asodesigner.fold import calculate_energies, get_weighted_energy
    from asodesigner.util import get_antisense

    all_data_human_gene_premrna_no_nan = df.copy()

    # Comment out the long cases for quick running
    for (window_size, step_size) in fold_variants:

        on_target_fold = 'on_target_fold_openness' + str(window_size) + '_' + str(step_size)
        on_target_fold_normalized = 'on_target_fold_openness_normalized' + str(window_size) + '_' + str(step_size)
        all_data_human_gene_premrna_no_nan[on_target_fold] = np.zeros_like(all_data_human_gene_premrna_no_nan[SEQUENCE],
                                                                           dtype=np.float64)
        all_data_human_gene_premrna_no_nan[on_target_fold_normalized] = np.zeros_like(
            all_data_human_gene_premrna_no_nan[SEQUENCE], dtype=np.float64)

        for gene in genes_u:

            target = gene_to_data[gene].full_mrna
            gene_rows = all_data_human_gene_premrna_no_nan[all_data_human_gene_premrna_no_nan[CANONICAL_GENE] == gene]
            energies = calculate_energies(str(target), step_size, window_size)

            for index, row in gene_rows.iterrows():
                antisense = row[SEQUENCE]
                sense = get_antisense(antisense)
                l = row[SENSE_LENGTH]
                sense_start = row[SENSE_START]
                mean_fold = get_weighted_energy(sense_start, l, step_size, energies, window_size)
                if mean_fold > 100:
                    print("Weird: ", mean_fold)
                    print("Sense_start ", sense_start)
                    print("Sense_length ", l)
                    print("Gene: ", gene)
                    mean_fold = 0
                all_data_human_gene_premrna_no_nan.loc[index, on_target_fold] = mean_fold
                all_data_human_gene_premrna_no_nan.loc[index, on_target_fold_normalized] = mean_fold / l
    return all_data_human_gene_premrna_no_nan


def populate_features(df, features):
    if 'self_energy' in features:
        df.loc[:, 'self_energy'] = [float(primer3.calc_homodimer(antisense).dg) for
                                    antisense in
                                    df[SEQUENCE]]
        df.loc[:, 'self_energy'] = df.loc[:,
                                   'self_energy'].astype(float)

    if 'internal_fold' in features:
        df.loc[:, 'internal_fold'] = [RNA.fold(antisense)[1] for antisense in df[SEQUENCE]]
