import json
import os
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict

import pandas as pd
from fuzzysearch import find_near_matches
from tqdm import tqdm

from asodesigner.cache import load_cache_off_target_hybridization, load_cache_off_target_wc, update_loaded_data, \
    save_cache
from asodesigner.consts import EXPERIMENT_RESULTS, CACHE_DIR
from asodesigner.experiment import Experiment
from asodesigner.features import SENSE_START, SENSE_LENGTH
from asodesigner.fold import get_trigger_mfe_scores_by_risearch, get_mfe_scores, dump_target_file, calculate_energies
from asodesigner.result import save_results_organism
from asodesigner.target_finder import iterate_template_antisense
from asodesigner.timer import Timer
from asodesigner.util import get_antisense


class LocusInfo:
    def __init__(self):
        self.exons = []
        self.introns = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.exon_concat = None
        self.full_mrna = None


def get_simplified_fasta_dict(fasta_dict):
    simplified_fasta_dict = dict()
    for locus_tag, locus_info in fasta_dict.items():
        simplified_fasta_dict[locus_tag] = str(locus_info.upper().seq)
    return simplified_fasta_dict


def validated_get_simplified_fasta_dict(fasta_dict, simplified_fasta_dict):
    if simplified_fasta_dict is None and fasta_dict is None:
        raise ValueError('Either simplified_fasta_dict or fasta_dict must be specified')

    if simplified_fasta_dict is None:
        return get_simplified_fasta_dict(fasta_dict)
    return simplified_fasta_dict


def process_fold_single_mrna(args):
    locus_tag, locus_info, step_size, window_size = args
    energies = calculate_energies(locus_info.full_mrna, step_size=step_size, window_size=window_size)
    return locus_tag, energies


def process_hybridization(task):
    i = task.sense_start
    l = task.sense_length
    locus_to_data = task.simplified_fasta_dict
    target_cache_filename = task.target_cache_filename

    parsing_type = task.parsing_type

    scores = get_trigger_mfe_scores_by_risearch(task.get_sense(), locus_to_data,
                                                minimum_score=task.minimum_score, neighborhood=l,
                                                parsing_type=parsing_type,
                                                target_file_cache=target_cache_filename)

    energy_scores = get_mfe_scores(scores, parsing_type=parsing_type)
    total_candidates = 0
    energy_sum = 0
    max_sum = 0
    binary_sum = 0
    for locus_scores in energy_scores:
        total_candidates += len(locus_scores)
        energy_sum += sum(locus_scores)

        if len(locus_scores) != 0:
            min_score = min(locus_scores)
        else:
            min_score = 0

        max_sum += min_score
        binary_sum += 1 if min_score < task.binary_cutoff else 0
    return ResultHybridization(i, l, total_candidates, energy_sum, max_sum, binary_sum)


def process_watson_crick_differences(args):
    i, l, aso_template, locus_to_data = args
    sense = aso_template[i:i + l]
    matches_per_distance = [0, 0, 0, 0]

    for locus_tag, locus_info in locus_to_data.items():
        matches = find_near_matches(sense, locus_info, max_insertions=0, max_deletions=0, max_l_dist=3)
        for match in matches:
            matches_per_distance[match.dist] += 1
            if match.dist == 0:
                print(locus_tag)

    # Return a tuple containing the starting index, current l, and match counts
    return (i, l, matches_per_distance[0],
            matches_per_distance[1], matches_per_distance[2], matches_per_distance[3])


def validate_organism(organism: str):
    organisms = ['human', 'yeast']
    if organism not in organisms:
        raise ValueError(f'Organism={organism} must be in {organisms}')


def parallelize_function(function, tasks, max_threads=None):
    """
    :param function: To be parallelized
    :param tasks: to be submitted to function
    :param max_threads: pass None to use all cores
    :return: results of parallel operation
    """
    results = []
    with Timer() as t:
        with ProcessPoolExecutor(max_workers=max_threads) as executor:
            futures = [executor.submit(function, task) for task in tasks]

            for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
                results.append(future.result())
    print(f"Parallel task done in: {t.elapsed_time}s")
    return results


def wc_results_to_dict(results, aso_template):
    results_dict = dict()
    for result in results:
        start = result[0]
        length = result[1]
        results_dict[get_antisense(aso_template[start:start + length])] = (result[2], result[3], result[4], result[5])

    return results_dict


def run_off_target_wc_analysis(experiment: Experiment, fasta_dict=None, simplified_fasta_dict=None, organism=None):
    validate_organism(organism)
    simplified_fasta_dict = validated_get_simplified_fasta_dict(fasta_dict, simplified_fasta_dict)

    loaded_data, cache_path = load_cache_off_target_wc(organism)
    aso_template = experiment.get_aso_template()

    tasks = []
    tasks_cached = 0
    for (i, l, antisense) in iterate_template_antisense(aso_template, experiment.l_values):
        if antisense not in loaded_data:  # no reason to calculate on cached elements
            tasks.append((i, l, experiment.get_aso_template(), simplified_fasta_dict))
        else:
            tasks_cached += 1

    print(f"Skipping {tasks_cached} tasks that were found in cache.")
    results = parallelize_function(process_watson_crick_differences, tasks)

    results_dict = wc_results_to_dict(results, aso_template)
    update_loaded_data(loaded_data, results_dict)
    save_cache(cache_path, loaded_data)

    full_results = []
    for i, l, antisense in iterate_template_antisense(aso_template, experiment.l_values):
        loaded_result = loaded_data[antisense]
        full_results.append((i, l, loaded_result[0], loaded_result[1], loaded_result[2], loaded_result[3]))

    columns = [SENSE_START, SENSE_LENGTH, '0_matches', '1_matches', '2_matches', '3_matches']
    df = pd.DataFrame(full_results, columns=columns)
    df = df.sort_values(by=['sense_start', 'sense_length'])

    print(df)
    save_results_organism(df, organism, experiment.name, 'wc_off_targets')


class Task:
    def __init__(self, sense_start, sense_length, aso_template, simplified_fasta_dict, target_cache_filename):
        self.sense_start = sense_start
        self.sense_length = sense_length
        self.aso_template = aso_template
        self.simplified_fasta_dict = simplified_fasta_dict
        self.target_cache_filename = target_cache_filename
        # Settings. TODO: consider moving to separate class
        self.minimum_score = 900
        self.parsing_type = '2'
        self.binary_cutoff = -20

    def get_sense(self):
        return self.aso_template[self.sense_start: self.sense_start + self.sense_length]

    def get_antisense(self):
        return get_antisense(self.get_sense())


@dataclass
class ResultHybridization:
    sense_start: int
    sense_length: int
    total_hybridization_candidates: int
    total_hybridization_energy: int
    total_hybridization_max_sum: int
    total_hybridization_binary_sum: int

    @staticmethod
    def results_to_result_dict(results, aso_template):
        results_dict = dict()
        for result in results:
            start = result.sense_start
            length = result.sense_length
            total_hybridization_candidates: int
            total_hybridization_energy: int
            total_hybridization_max_sum: int
            total_hybridization_binary_sum: int

            results_dict[get_antisense(aso_template[start:start + length])] = (
                result.total_hybridization_candidates, result.total_hybridization_energy,
                result.total_hybridization_max_sum, result.total_hybridization_binary_sum)

        return results_dict


def run_off_target_hybridization_analysis(experiment: Experiment, fasta_dict=None, simplified_fasta_dict=None,
                                          organism=None):
    validate_organism(organism)
    simplified_fasta_dict = validated_get_simplified_fasta_dict(fasta_dict, simplified_fasta_dict)

    hash = random.getrandbits(64)

    target_cache_filename = f'target-cache-{hash}.fa'
    # to improve speed of process_hybridization
    target_cache_path = dump_target_file(target_cache_filename, simplified_fasta_dict)

    loaded_data, cache_path = load_cache_off_target_hybridization(organism)
    aso_template = experiment.get_aso_template()

    tasks = []
    tasks_cached = 0

    for i, l, antisense in iterate_template_antisense(aso_template, experiment.l_values):
        if not antisense in loaded_data:
            tasks.append(Task(i, l, aso_template, simplified_fasta_dict, str(target_cache_path)))
        else:
            tasks_cached += 1
    print(f"Skipping {tasks_cached} tasks that were found in cache.")

    results = parallelize_function(process_hybridization, tasks)

    results_dict = ResultHybridization.results_to_result_dict(results, aso_template)

    update_loaded_data(loaded_data, results_dict)
    save_cache(cache_path, loaded_data)

    full_results = []
    for i, l, antisense in iterate_template_antisense(aso_template, experiment.l_values):
        loaded_result = loaded_data[antisense]
        full_results.append((i, l, loaded_result[0], loaded_result[1], loaded_result[2], loaded_result[3]))

    columns = ['sense_start', 'sense_length', 'total_hybridization_candidates', 'total_hybridization_energy',
               'total_hybridization_max_sum', 'total_hybridization_binary_sum']
    df = pd.DataFrame([result for result in full_results], columns=columns)

    print(df)
    df = df.sort_values(by=['sense_start', 'sense_length'])
    save_results_organism(df, organism, experiment.name, 'hybridization_off_targets')
    os.remove(target_cache_path)


def run_off_target_fold_analysis(locus_to_data, experiment_name, organism):
    validate_organism(organism)

    window_size = 40
    step_size = 15

    tasks = []
    for key, value in locus_to_data.items():
        tasks.append((key, value, step_size, window_size))

    results = parallelize_function(process_fold_single_mrna, tasks)
    results_dict = dict()
    for result in results:
        results_dict[result[0]] = result[1]

    result_path = (EXPERIMENT_RESULTS / experiment_name /
                   f"{organism}_results" / f'fold_off_target_fold_energy_window_{window_size}_step_{step_size}.json')
    with open(result_path, 'w') as f:
        json.dump(results, f, indent=4)
