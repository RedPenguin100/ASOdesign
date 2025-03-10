import pandas as pd

from Bio import SeqIO
from fuzzysearch import find_near_matches

from asodesigner.util import get_antisense


def get_gfp_seq_and_context():
    gfp_obj = next(SeqIO.parse('../data/GFP_first_exp.fasta', 'fasta'))
    gfp_seq = str(gfp_obj.seq.upper())

    with open('../data/GFP_context.txt', 'r') as f:
        gfp_context = f.read().upper()
    return (gfp_seq, gfp_context)

def get_gfp_first_exp(gap=100):
    # TODO: gap should be always 100 in this function
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    gfp_ext = gfp_context[gfp_start - gap: gfp_start + len(gfp_seq) + gap]

    return gfp_ext


def get_gfp_second_exp():
    right_gap = 50
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    gfp_ext = gfp_context[gfp_start: gfp_start + len(gfp_seq) + right_gap]
    return gfp_ext

def generate_scrambled(target_seq):
    l_values = [17, 18, 19, 20, 21]

    matches_per_distance = [0, 0, 0, 0]
    df = pd.DataFrame(
        columns=['sense_start', 'sense_length', '0_matches', '1_matches', '2_matches', '3_matches'])

    for l in l_values:
        for i in range(len(target_seq) - l + 1):
            sense = target_seq[i:i + l]
            sample = get_antisense(sense)

            matches = find_near_matches(sample, target_seq, max_insertions=0, max_deletions=0, max_l_dist=3)
            for match in matches:
                matches_per_distance[match.dist] += 1

            if (i, l) not in df.index:
                df.loc[len(df)] = [i, l, matches_per_distance[0], matches_per_distance[1], matches_per_distance[2],
                                   matches_per_distance[3]]
    return df


if __name__ == '__main__':
    df = generate_scrambled(get_gfp_first_exp())
    print(df[df['0_matches'] != 0])
    print(df[df['1_matches'] != 0])
    print(df[df['2_matches'] != 0])
    print(df[df['3_matches'] != 0])

