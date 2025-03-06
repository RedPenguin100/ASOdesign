import math

from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from fuzzysearch import find_near_matches

from asodesigner.consts import HUMAN_TRANSCRIPTS_FASTA
from asodesigner.target_finder import get_gfp_first_exp
from asodesigner.timer import Timer
from asodesigner.util import get_antisense


def process_sense(args):
    i, l, gfp_seq, fasta_dict = args
    sense = gfp_seq[i:i + l]
    matches_per_distance = [0, 0, 0, 0]

    for locus_tag, locus_info in fasta_dict.items():
        off_target = str(locus_info.upper().seq)
        matches = find_near_matches(sense, off_target, max_insertions=0, max_deletions=0, max_l_dist=3)
        for match in matches:
            matches_per_distance[match.dist] += 1
            if match.dist == 0:
                print(locus_tag)

    # Return a tuple containing the starting index, current l, and match counts
    return (i, l, matches_per_distance[0],
            matches_per_distance[1], matches_per_distance[2], matches_per_distance[3])


class Experiment:
    def __init__(self):
        self.name = None
        self.target_sequence = None


def main():
    with Timer() as t:
        with open(str(HUMAN_TRANSCRIPTS_FASTA), 'r') as handle:
            fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    print(f"Reading FASTA took: {t.elapsed_time}s")

    tasks = []
    # gfp_seq = get_antisense(get_gfp_first_exp())

    experiment = Experiment()
    experiment.name = "FirstScrambled"
    experiment.target_sequence = get_antisense(get_gfp_first_exp())

    print(experiment.target_sequence)
    l_values = [16, 17, 18, 19, 20, 21, 22]
    max_iterations = math.inf

    for l in l_values:
        for i in range(min(len(experiment.target_sequence) - l + 1, max_iterations)):
            tasks.append((i, l, experiment.target_sequence, fasta_dict))

    results = []
    with Timer() as t:
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_sense, arg) for arg in tasks]

            for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
                results.append(future.result())
    print(f"Read off targets in: {t.elapsed_time}s")

    columns = ['sense_start', 'sense_length', '0_matches', '1_matches', '2_matches', '3_matches']
    df = pd.DataFrame(results, columns=columns)

    print(df)
    df.to_csv(f'human_results/{experiment.name}gfp_off_targets.csv', index=False)


if __name__ == "__main__":
    main()
