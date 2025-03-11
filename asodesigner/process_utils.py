import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from fuzzysearch import find_near_matches
from tqdm import tqdm

from asodesigner.experiment import Experiment
from asodesigner.fold import get_trigger_mfe_scores_by_risearch, get_mfe_scores, dump_target_file
from asodesigner.timer import Timer


def get_simplified_fasta_dict(fasta_dict):
    simplified_fasta_dict = dict()
    for locus_tag, locus_info in fasta_dict.items():
        simplified_fasta_dict[locus_tag] = str(locus_info.upper().seq)
    return simplified_fasta_dict


def process_hybridization(args):
    i, l, target_seq, locus_to_data, target_cache_filename = args

    parsing_type = '2'
    scores = get_trigger_mfe_scores_by_risearch(target_seq[i: i + l], locus_to_data,
                                                minimum_score=900, neighborhood=l, parsing_type=parsing_type,
                                                target_file_cache=target_cache_filename)
    energy_scores = get_mfe_scores(scores, parsing_type=parsing_type)
    total_candidates = 0
    energy_sum = 0
    max_sum = 0
    binary_sum = 0
    for locus_scores in energy_scores:
        total_candidates += len(locus_scores)
        energy_sum += sum(locus_scores)
        max_sum += min(locus_scores)
        binary_sum += 1 if min(locus_scores) < -20 else 0
    return (i, l, total_candidates, energy_sum, max_sum, binary_sum)


def process_watson_crick_differences(args):
    i, l, target_seq, locus_to_data = args
    sense = target_seq[i:i + l]
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


def run_off_target_wc_analysis(experiment: Experiment, fasta_dict=None, simplified_fasta_dict=None, organism=None):
    if organism not in ['human', 'yeast']:
        raise ValueError('Organism must be either "human" or "yeast"')

    if simplified_fasta_dict is None and fasta_dict is None:
        raise ValueError('Either simplified_fasta_dict or fasta_dict must be specified')

    if simplified_fasta_dict is None:
        simplified_fasta_dict = get_simplified_fasta_dict(fasta_dict)

    tasks = []
    for l in experiment.l_values:
        for i in range(len(experiment.target_sequence) - l + 1):
            tasks.append((i, l, experiment.target_sequence, simplified_fasta_dict))

    results = []
    with Timer() as t:
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_watson_crick_differences, arg) for arg in tasks]

            for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
                results.append(future.result())
    print(f"Read off targets in: {t.elapsed_time}s")

    columns = ['sense_start', 'sense_length', '0_matches', '1_matches', '2_matches', '3_matches']
    df = pd.DataFrame(results, columns=columns)

    print(df)
    df.to_csv(f'{organism}_results/{experiment.name}gfp_off_targets.csv', index=False)



def run_off_target_hybridization_analysis(experiment: Experiment, fasta_dict=None, simplified_fasta_dict=None,
                                          organism=None):
    if organism not in ['human', 'yeast']:
        raise ValueError('Organism must be either "human" or "yeast"')

    if fasta_dict is None and simplified_fasta_dict is None:
        raise ValueError('Either simplified_fasta_dict or fasta_dict must be specified')

    if simplified_fasta_dict is None:
        simplified_fasta_dict = get_simplified_fasta_dict(fasta_dict)

    target_cache_filename = 'target-cache.fa'

    tasks = []
    for l in experiment.l_values:
        for i in range(len(experiment.target_sequence) - l + 1):
            tasks.append((i, l, experiment.target_sequence, simplified_fasta_dict, target_cache_filename))

    dump_target_file(target_cache_filename, simplified_fasta_dict)

    results = []
    with Timer() as t:
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_hybridization, arg) for arg in tasks]

            for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
                results.append(future.result())
    print(f"Read Hybridization off targets in: {t.elapsed_time}s")

    columns = ['sense_start', 'sense_length', 'total_hybridization_candidates', 'total_hybridization_energy',
               'total_hybridization_max_sum', 'total_hybridization_binary_sum']
    df = pd.DataFrame(results, columns=columns)

    print(df)
    df.to_csv(f'{organism}_results/{experiment.name}gfp_hybridization_off_targets.csv', index=False)

    os.remove(target_cache_filename)