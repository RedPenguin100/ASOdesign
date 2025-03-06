import time
import math

from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from fuzzysearch import find_near_matches

from asodesigner.consts import HUMAN_TRANSCRIPTS_FASTA
from asodesigner.target_finder import get_gfp_first_exp


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


def main():
    start = time.time()
    with open(str(HUMAN_TRANSCRIPTS_FASTA), 'r') as handle:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    end = time.time()
    print(f"Reading FASTA took: {end - start}s")

    tasks = []
    gfp_seq = get_gfp_first_exp()
    print(gfp_seq)
    # l_values = [15, 16, 17, 18, 19, 20, 21, 22]
    # max_iterations = math.inf
    #
    # for l in l_values:
    #     for i in range(min(len(gfp_seq) - l + 1, max_iterations)):
    #         tasks.append((i, l, gfp_seq, fasta_dict))

    l_values = [19]

    max_iterations = math.inf

    for i in range(50, 100):
        tasks.append((i, 19, gfp_seq, fasta_dict))

    results = []
    start = time.time()
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_sense, arg) for arg in tasks]

        for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
            results.append(future.result())
    end = time.time()
    print(f"Read off targets in: {end - start}s")

    columns = ['sense_start', 'sense_length', '0_matches', '1_matches', '2_matches', '3_matches']
    df = pd.DataFrame(results, columns=columns)

    print(f"Time took to find: {end - start}s")

    print(df)
    df.to_csv('human_results/gfp_off_targets.csv', index=False)


if __name__ == "__main__":
    main()
