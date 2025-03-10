import math
import os

from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

from asodesigner.consts import HUMAN_TRANSCRIPTS_FASTA
from asodesigner.fold import dump_target_file
from asodesigner.read_utils import process_hybridization, process_watson_crick_differences
from asodesigner.target_finder import get_gfp_second_exp
from asodesigner.timer import Timer
from asodesigner.util import get_antisense


class Experiment:
    def __init__(self):
        self.name = None
        self.target_sequence = None


def get_simplified_fasta_dict(fasta_dict):
    simplified_fasta_dict = dict()
    for locus_tag, locus_info in fasta_dict.items():
        simplified_fasta_dict[locus_tag] = str(locus_info.upper().seq)
    return simplified_fasta_dict


def run_off_target_wc_analysis(fasta_dict, experiment):
    l_values = [16, 17, 18, 19, 20, 21, 22]

    simplified_fasta_dict = get_simplified_fasta_dict(fasta_dict)

    tasks = []
    for l in l_values:
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
    df.to_csv(f'human_results/{experiment.name}gfp_off_targets.csv', index=False)


def run_off_target_hybridization_analysis(fasta_dict, experiment):
    l_values = [16, 17, 18, 19, 20, 21, 22]

    simplified_fasta_dict = get_simplified_fasta_dict(fasta_dict)

    target_cache_filename = 'target-cache.fa'

    tasks = []
    for l in l_values:
        for i in range(len(experiment.target_sequence) - l + 1):
            tasks.append((i, l, experiment.target_sequence, simplified_fasta_dict, target_cache_filename))

    dump_target_file(target_cache_filename, simplified_fasta_dict)

    results = []
    with Timer() as t:
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_hybridization, arg) for arg in tasks]

            for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
                results.append(future.result())
    print(f"Read off targets in: {t.elapsed_time}s")

    columns = ['sense_start', 'sense_length', 'total_hybridization_candidates', 'total_hybridization_energy',
               'total_hybridization_max_sum', 'total_hybridization_binary_sum']
    df = pd.DataFrame(results, columns=columns)

    print(df)
    df.to_csv(f'human_results/{experiment.name}gfp_hybridization_off_targets.csv', index=False)

    os.remove(target_cache_filename)


def main():
    with Timer() as t:
        with open(str(HUMAN_TRANSCRIPTS_FASTA), 'r') as handle:
            fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    print(f"Reading FASTA took: {t.elapsed_time}s")

    experiment = Experiment()
    experiment.name = "SecondScrambled"
    experiment.target_sequence = get_antisense(get_gfp_second_exp())

    print(experiment.target_sequence)

    run_off_target_hybridization_analysis(fasta_dict, experiment)


if __name__ == "__main__":
    main()
