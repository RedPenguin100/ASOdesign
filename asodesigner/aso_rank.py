import numpy as np
import pandas as pd
from ViennaRNA import RNA
import math
from Bio.SeqUtils import gc_fraction, MeltingTemp

from asodesigner.fold import get_weighted_energy, calculate_energies
from asodesigner.target_finder import get_gfp_first_exp
from asodesigner.util import get_antisense


def record_internal_fold(target_seq, l_values, experiment_name):
    max_iterations = math.inf
    results = []
    for l in l_values:
        for i in range(min(len(target_seq) - l + 1, max_iterations)):
            sense = target_seq[i:i + l]
            antisense = get_antisense(sense)
            structure, mfe = RNA.fold(antisense)
            results.append((i, l, structure, mfe))

    columns = ['sense_start', 'sense_length', 'structure', 'mfe']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv('antisense_results/' + experiment_name + 'antisense_fold.csv', index=False)
    return df


def record_nucleotide_properties(target_seq, l_values, experiment_name):
    max_iterations = math.inf
    results = []
    for l in l_values:
        for i in range(min(len(target_seq) - l + 1, max_iterations)):
            sense = target_seq[i:i + l]
            antisense = get_antisense(sense)
            gc_content = gc_fraction(antisense)
            contains_GGGG = 'GGGG' in antisense
            results.append((i, l, gc_content, contains_GGGG))
    columns = ['sense_start', 'sense_length', 'gc_content', 'contains_GGGG']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv('antisense_results/' + experiment_name + 'antisense_nucleotide_properties.csv', index=False)
    return df


def record_melting_temperature(target_seq, l_values, experiment_name):
    max_iterations = math.inf
    results = []
    for l in l_values:
        for i in range(min(len(target_seq) - l + 1, max_iterations)):
            sense = target_seq[i:i + l]
            antisense = get_antisense(sense)
            melting_temperature = MeltingTemp.Tm_NN(antisense)

            results.append((i, l, melting_temperature))
    columns = ['sense_start', 'sense_length', 'melting_temperature']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv('antisense_results/' + experiment_name + 'antisense_melting_temperature.csv', index=False)
    return df




def record_on_target_fold(target_seq, l_values, experiment_name):
    window_size = 40
    step_size = 15

    energies = calculate_energies(target_seq, step_size, window_size)

    results = []
    for l in l_values:
        for i in range(len(target_seq) - l + 1):
            mean_fold = get_weighted_energy(i, l, step_size, energies, window_size)
            results.append((i, l, mean_fold, mean_fold / l))

    columns = ['sense_start', 'sense_length', 'fold_openness', 'fold_openness_normalized']

    df = pd.DataFrame(results, columns=columns)
    return df


if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 1000)

    l_values = [15, 16, 17, 18, 19, 20, 21]

    gfp_seq = get_gfp_first_exp(gap=0)
    print("Target length: ", len(gfp_seq))
    experiment_name = 'Trying'

    # record_internal_fold(gfp_seq, l_values, experiment_name)
    # record_nucleotide_properties(gfp_seq, l_values, experiment_name)
    # df = record_melting_temperature(gfp_seq, l_values, experiment_name)
    df = record_on_target_fold(gfp_seq, l_values, experiment_name)
    print(df)

    # print(df.sort_values(by='melting_temperature'))
