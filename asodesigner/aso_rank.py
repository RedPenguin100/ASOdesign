import math
import pandas as pd
from ViennaRNA import RNA
import primer3

from Bio.SeqUtils import gc_fraction, MeltingTemp
from fuzzysearch import find_near_matches

from asodesigner.experiment import Experiment, get_experiments, maybe_create_experiment_folders
from asodesigner.fold import get_weighted_energy, calculate_energies, get_trigger_mfe_scores_by_risearch, get_mfe_scores
from asodesigner.result import save_results_on_target
from asodesigner.timer import Timer
from asodesigner.util import get_antisense


def record_internal_fold(experiment: Experiment):
    results = []
    for (i, l, antisense) in experiment.get_aso_antisense_iterator():
        structure, mfe = RNA.fold(antisense)
        results.append((i, l, structure, mfe))

    columns = ['sense_start', 'sense_length', 'structure', 'mfe']

    df = pd.DataFrame(results, columns=columns)
    save_results_on_target(df, experiment.name, 'antisense_fold')
    return df


def record_nucleotide_properties(experiment: Experiment):
    results = []
    for (i, l, antisense) in experiment.get_aso_antisense_iterator():
        gc_content = gc_fraction(antisense)
        contains_GGGG = 'GGGG' in antisense
        results.append((i, l, gc_content, contains_GGGG))
    columns = ['sense_start', 'sense_length', 'gc_content', 'contains_GGGG']

    df = pd.DataFrame(results, columns=columns)
    save_results_on_target(df, experiment.name, 'antisense_nucleotide_properties')
    return df


def record_melting_temperature(experiment: Experiment):
    results = []
    for (i, l, antisense) in experiment.get_aso_antisense_iterator():
        melting_temperature = MeltingTemp.Tm_NN(antisense)

        results.append((i, l, melting_temperature))
    columns = ['sense_start', 'sense_length', 'melting_temperature']

    df = pd.DataFrame(results, columns=columns)
    save_results_on_target(df, experiment.name, 'antisense_melting_temperature')
    return df


def record_self_dimerization_unmodified(experiment: Experiment):
    results = []
    for (i, l, antisense) in experiment.get_aso_antisense_iterator():
        homodimer = primer3.calc_homodimer(antisense)
        results.append((i, l, homodimer.dg, homodimer.dh, homodimer.ds, homodimer.structure_found, homodimer.tm))

    columns = ['sense_start', 'sense_length', 'aso_aso_delta_g', 'aso_aso_delta_h', 'aso_aso_delta_s',
               'aso_aso_structure_found', 'aso_aso_tm']

    df = pd.DataFrame(results, columns=columns)
    save_results_on_target(df, experiment.name, 'antisense_self_dimerization_unmodified')
    return df


def record_on_target_fold(experiment: Experiment):
    window_size = 40
    step_size = 15

    energies = calculate_energies(experiment.target_sequence, step_size, window_size)

    results = []
    for (i, l, _) in experiment.get_aso_sense_iterator():
        mean_fold = get_weighted_energy(i, l, step_size, energies, window_size)
        results.append((i, l, mean_fold, mean_fold / l))

    columns = ['sense_start', 'sense_length', 'on_target_fold_openness', 'on_target_fold_openness_normalized']

    df = pd.DataFrame(results, columns=columns)
    save_results_on_target(df, experiment.name, 'on_target_fold')
    return df


def record_on_target_energy_hybridization(experiment: Experiment):
    name_to_sequence = {'target_seq': experiment.target_sequence}
    results = []
    parsing_type = '2'
    for (i, l, sense) in experiment.get_aso_sense_iterator():
        tmp_results = get_trigger_mfe_scores_by_risearch(sense, name_to_sequence, minimum_score=1200, neighborhood=l,
                                                         parsing_type=parsing_type)
        scores = get_mfe_scores(tmp_results, parsing_type)
        if len(scores) == 0:
            results.append((i, l, 0, 0., 0.))
        else:
            target_scores = scores[0]
            min_score = 0 if len(target_scores) == 0 else min(target_scores)
            results.append((i, l, len(target_scores), sum(target_scores), min_score))

    column = ['sense_start', 'sense_length', 'on_target_energy_fits', 'on_target_energy_sum', 'on_target_energy_max']
    df = pd.DataFrame(results, columns=column)
    save_results_on_target(df, experiment.name, 'on_target_energy')
    return df


def record_on_target_wc_hybridization(experiment: Experiment):
    results = []

    max_distance = 3

    for (idx, l, antisense) in experiment.get_aso_antisense_iterator():
        matches_per_distance = [0 for it in range(max_distance + 1)]
        matches = find_near_matches(antisense, experiment.target_sequence, max_insertions=0, max_deletions=0,
                                    max_l_dist=max_distance)
        for match in matches:
            matches_per_distance[match.dist] += 1

        result = [idx, l]
        for distance in range(max_distance + 1):
            result.append(matches_per_distance[distance])

        results.append(tuple(result))

    columns = ['sense_start', 'sense_length']
    for it in range(max_distance + 1):
        columns.append(f'on_matches_{it}')
    df = pd.DataFrame(results, columns=columns)
    save_results_on_target(df, experiment.name, 'on_target_wc')
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


    experiment_names = ['EntirePositiveControl']
    experiments = get_experiments(experiment_names)

    for experiment in experiments:
        maybe_create_experiment_folders(experiment.name)

        print("Target length: ", len(experiment.target_sequence))

        with Timer() as t:
            all_record(experiment)

        print(f"Recording all {experiment.name} took: ", t.elapsed_time)
