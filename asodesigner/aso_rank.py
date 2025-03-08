import numpy as np
import pandas as pd
from ViennaRNA import RNA
import math
import primer3


from Bio.SeqUtils import gc_fraction, MeltingTemp
from fuzzysearch import find_near_matches

from asodesigner.fold import get_weighted_energy, calculate_energies, get_trigger_mfe_scores_by_risearch, get_mfe_scores
from asodesigner.target_finder import get_gfp_first_exp
from asodesigner.util import get_antisense
from numba import njit


@njit
def iterate_target(target_seq, l_values):
    for l in l_values:
        for i in range(len(target_seq) - l + 1):
            yield i, l, target_seq[i:i + l]


def record_internal_fold(target_seq, l_values, experiment_name):
    results = []
    for (i, l, sense) in iterate_target(target_seq, l_values):
        antisense = get_antisense(sense)
        structure, mfe = RNA.fold(antisense)
        results.append((i, l, structure, mfe))

    columns = ['sense_start', 'sense_length', 'structure', 'mfe']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment_name}_antisense_fold.csv', index=False)
    return df


def record_nucleotide_properties(target_seq, l_values, experiment_name):
    results = []
    for (i, l, sense) in iterate_target(target_seq, l_values):
        antisense = get_antisense(sense)
        gc_content = gc_fraction(antisense)
        contains_GGGG = 'GGGG' in antisense
        results.append((i, l, gc_content, contains_GGGG))
    columns = ['sense_start', 'sense_length', 'gc_content', 'contains_GGGG']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment_name}antisense_nucleotide_properties.csv', index=False)
    return df


def record_melting_temperature(target_seq, l_values, experiment_name):
    results = []
    for (i, l, sense) in iterate_target(target_seq, l_values):
        antisense = get_antisense(sense)
        melting_temperature = MeltingTemp.Tm_NN(antisense)

        results.append((i, l, melting_temperature))
    columns = ['sense_start', 'sense_length', 'melting_temperature']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment_name}antisense_melting_temperature.csv', index=False)
    return df


def record_self_dimerization_unmodified(target_seq, l_values, experiment_name):
    results = []
    for (i, l, sense) in iterate_target(target_seq, l_values):
        antisense = get_antisense(sense)
        homodimer = primer3.calc_homodimer(antisense)
        results.append((i, l, homodimer.dg, homodimer.dh, homodimer.ds, homodimer.structure_found, homodimer.tm))

    columns = ['sense_start', 'sense_length', 'delta_g', 'delta_h', 'delta_s', 'structure_found', 'tm']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment_name}_antisense_self_dimerization_unmodified.csv', index=False)
    return df


def record_on_target_fold(target_seq, l_values, experiment_name):
    window_size = 40
    step_size = 15

    energies = calculate_energies(target_seq, step_size, window_size)

    results = []
    for (i, l, sense) in iterate_target(target_seq, l_values):
        mean_fold = get_weighted_energy(i, l, step_size, energies, window_size)
        results.append((i, l, mean_fold, mean_fold / l))

    columns = ['sense_start', 'sense_length', 'fold_openness', 'fold_openness_normalized']

    df = pd.DataFrame(results, columns=columns)
    df.to_csv(f'antisense_results/{experiment_name}on_target_fold.csv', index=False)

    return df


def record_on_target_energy_hybridization(target_seq, l_values, experiment_name):
    name_to_sequence = {'target_seq': target_seq}
    results = []
    for (i, l, sense) in iterate_target(target_seq, l_values):
        tmp_results = get_trigger_mfe_scores_by_risearch(sense, name_to_sequence)
        scores = get_mfe_scores(tmp_results)
        print(scores)
        results.append((i, l, len(scores[0])))

    column = ['sense_start', 'sense_length', 'energy_fits']
    df = pd.DataFrame(results, columns=column)
    return df


def record_on_target_wc_hybridization(target_seq, l_values, experiment_name):
    results = []

    max_distance = 3

    for (i, l, sense) in iterate_target(target_seq, l_values):
        matches_per_distance = [0 for i in range(max_distance + 1)]
        matches = find_near_matches(get_antisense(sense), target_seq, max_insertions=0, max_deletions=0,
                                    max_l_dist=max_distance)
        for match in matches:
            matches_per_distance[match.dist] += 1

        result = [i, l]
        for distance in range(max_distance + 1):
            result.append(matches_per_distance[distance])

        results.append(tuple(result))

    columns = ['sense_start', 'sense_length']
    for i in range(max_distance + 1):
        columns.append(f'on_matches_{i}')
    df = pd.DataFrame(results, columns=columns)
    return df


if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 1000)

    l_values = [15, 16, 17, 18, 19, 20, 21]
    # l_values = [15]

    gfp_seq = get_gfp_first_exp(gap=100)
    print(gfp_seq[:15])
    print(gfp_seq[1:16])
    print(gfp_seq[3:18])
    print(gfp_seq[81:101])
    print(gfp_seq[127:127+15])

    print("Target length: ", len(gfp_seq))
    experiment_name = 'First'

    # record_internal_fold(gfp_seq, l_values, experiment_name)
    # record_on_target_fold(gfp_seq, l_values, experiment_name)
    record_self_dimerization_unmodified(gfp_seq, l_values, experiment_name)
    # record_nucleotide_properties(gfp_seq, l_values, experiment_name)
    # df = record_melting_temperature(gfp_seq, l_values, experiment_name)
    # df = record_on_target_fold(gfp_seq, l_values, experiment_name)
    # df = record_on_target_energy_hybridization(gfp_seq, l_values, experiment_name)
    # print(df)
    # df = record_on_target_wc_hybridization(gfp_seq, l_values, experiment_name)
    # print(df)
    # print(df.sort_values(by='melting_temperature'))
