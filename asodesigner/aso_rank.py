import numpy as np
import pandas as pd
from ViennaRNA import RNA
import math
import primer3

from Bio.SeqUtils import gc_fraction, MeltingTemp
from fuzzysearch import find_near_matches

from asodesigner.experiment import Experiment
from asodesigner.fold import get_weighted_energy, calculate_energies, get_trigger_mfe_scores_by_risearch, get_mfe_scores
from asodesigner.target_finder import get_gfp_first_exp, get_gfp_second_exp
from asodesigner.timer import Timer
from asodesigner.util import get_antisense

from numba import njit


@njit
def iterate_target(target_seq, l_values):
    for l in l_values:
        for i in range(len(target_seq) - l + 1):
            yield i, l, target_seq[i:i + l]


def record_internal_fold(experiment: Experiment):
    results = []
    for (i, l, sense) in iterate_target(experiment.target_sequence, experiment.l_values):
        antisense = get_antisense(sense)
        structure, mfe = RNA.fold(antisense)
        results.append((i, l, structure, mfe))

    columns = ['sense_start', 'sense_length', 'structure', 'mfe']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment.name}_antisense_fold.csv', index=False)
    return df


def record_nucleotide_properties(experiment: Experiment):
    results = []
    for (i, l, sense) in iterate_target(experiment.target_sequence, experiment.l_values):
        antisense = get_antisense(sense)
        gc_content = gc_fraction(antisense)
        contains_GGGG = 'GGGG' in antisense
        results.append((i, l, gc_content, contains_GGGG))
    columns = ['sense_start', 'sense_length', 'gc_content', 'contains_GGGG']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment.name}antisense_nucleotide_properties.csv', index=False)
    return df


def record_melting_temperature(experiment: Experiment):
    results = []
    for (i, l, sense) in iterate_target(experiment.target_sequence, experiment.l_values):
        antisense = get_antisense(sense)
        melting_temperature = MeltingTemp.Tm_NN(antisense)

        results.append((i, l, melting_temperature))
    columns = ['sense_start', 'sense_length', 'melting_temperature']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment.name}antisense_melting_temperature.csv', index=False)
    return df


def record_self_dimerization_unmodified(experiment: Experiment):
    results = []
    for (i, l, sense) in iterate_target(experiment.target_sequence, experiment.l_values):
        antisense = get_antisense(sense)
        homodimer = primer3.calc_homodimer(antisense)
        results.append((i, l, homodimer.dg, homodimer.dh, homodimer.ds, homodimer.structure_found, homodimer.tm))

    columns = ['sense_start', 'sense_length', 'delta_g', 'delta_h', 'delta_s', 'structure_found', 'tm']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv(f'antisense_results/{experiment.name}_antisense_self_dimerization_unmodified.csv', index=False)
    return df


def record_on_target_fold(experiment: Experiment):
    window_size = 40
    step_size = 15

    energies = calculate_energies(experiment.target_sequence, step_size, window_size)

    results = []
    for (i, l, sense) in iterate_target(experiment.target_sequence, experiment.l_values):
        mean_fold = get_weighted_energy(i, l, step_size, energies, window_size)
        results.append((i, l, mean_fold, mean_fold / l))

    columns = ['sense_start', 'sense_length', 'fold_openness', 'fold_openness_normalized']

    df = pd.DataFrame(results, columns=columns)
    df.to_csv(f'antisense_results/{experiment.name}on_target_fold.csv', index=False)

    return df


def record_on_target_energy_hybridization(experiment: Experiment):
    name_to_sequence = {'target_seq': experiment.target_sequence}
    results = []
    parsing_type = '2'
    for (i, l, sense) in iterate_target(experiment.target_sequence, experiment.l_values):
        tmp_results = get_trigger_mfe_scores_by_risearch(sense, name_to_sequence, minimum_score=900, neighborhood=l,
                                                         parsing_type=parsing_type)
        scores = get_mfe_scores(tmp_results, parsing_type)
        if len(scores) == 0:
            results.append((i, l, 0, 0., 0.))
        else:
            target_scores = scores[0]
            results.append((i, l, len(target_scores), sum(target_scores), max(target_scores)))

    column = ['sense_start', 'sense_length', 'on_target_energy_fits', 'on_target_energy_sum', 'on_target_energy_max']
    df = pd.DataFrame(results, columns=column)
    df.to_csv(f'antisense_results/{experiment.name}on_target_energy.csv', index=False)
    return df


def record_on_target_wc_hybridization(experiment: Experiment):
    results = []

    max_distance = 3

    for (i, l, sense) in iterate_target(experiment.target_sequence, experiment.l_values):
        matches_per_distance = [0 for i in range(max_distance + 1)]
        matches = find_near_matches(get_antisense(sense), experiment.target_sequence, max_insertions=0, max_deletions=0,
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
    df.to_csv(f'antisense_results/{experiment.name}on_target_wc.csv', index=False)
    return df


def all_record(experiment: Experiment):
    # Fold properties
    record_internal_fold(experiment)
    record_on_target_fold(experiment)
    record_self_dimerization_unmodified(experiment)

    # Energy properties
    record_melting_temperature(experiment)
    record_on_target_energy_hybridization(experiment)

    # Nucleotide properties
    record_nucleotide_properties(experiment)
    # record_on_target_wc_hybridization(experiment)


if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 1000)

    experiment = Experiment()
    experiment.target_sequence = get_antisense(get_gfp_second_exp())
    experiment.name = 'SecondScrambled'
    experiment.l_values = [16, 17, 18, 19, 20, 21, 22]

    print("Target length: ", len(experiment.target_sequence))

    with Timer() as t:
        all_record(experiment)

    print("Recording all took: ", t.elapsed_time)
