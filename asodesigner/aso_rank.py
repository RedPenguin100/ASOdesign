import pandas as pd
from ViennaRNA import RNA
import math
from Bio.SeqUtils import gc_fraction, MeltingTemp

from asodesigner.target_finder import get_gfp_first_exp
from asodesigner.util import get_antisense


def record_internal_fold(gfp_seq, l_values, experiment_name):
    max_iterations = math.inf
    results = []
    for l in l_values:
        for i in range(min(len(gfp_seq) - l + 1, max_iterations)):
            sense = gfp_seq[i:i + l]
            antisense = get_antisense(sense)
            structure, mfe = RNA.fold(antisense)
            results.append((i, l, structure, mfe))

    columns = ['sense_start', 'sense_length', 'structure', 'mfe']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv('antisense_results/' + experiment_name + 'antisense_fold.csv', index=False)
    return df


def record_nucleotide_properties(gfp_seq, l_values, experiment_name):
    max_iterations = math.inf
    results = []
    for l in l_values:
        for i in range(min(len(gfp_seq) - l + 1, max_iterations)):
            sense = gfp_seq[i:i + l]
            antisense = get_antisense(sense)
            gc_content = gc_fraction(antisense)
            contains_GGGG = 'GGGG' in antisense
            results.append((i, l, gc_content, contains_GGGG))
    columns = ['sense_start', 'sense_length', 'gc_content', 'contains_GGGG']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv('antisense_results/' + experiment_name + 'antisense_nucleotide_properties.csv', index=False)
    return df


def record_melting_temperature(gfp_seq, l_values, experiment_name):
    max_iterations = math.inf
    results = []
    for l in l_values:
        for i in range(min(len(gfp_seq) - l + 1, max_iterations)):
            sense = gfp_seq[i:i + l]
            antisense = get_antisense(sense)
            melting_temperature = MeltingTemp.Tm_NN(antisense)

            results.append((i, l, melting_temperature))
    columns = ['sense_start', 'sense_length', 'melting_temperature']

    df = pd.DataFrame(results, columns=columns)

    df.to_csv('antisense_results/' + experiment_name + 'antisense_melting_temperature.csv', index=False)
    return df


if __name__ == '__main__':
    l_values = [15, 16, 17, 18, 19, 20, 21]

    gfp_seq = get_gfp_first_exp()
    experiment_name = 'First'

    record_internal_fold(gfp_seq, l_values, experiment_name)
    record_nucleotide_properties(gfp_seq, l_values, experiment_name)
    df = record_melting_temperature(gfp_seq, l_values, experiment_name)
    print(df.sort_values(by='melting_temperature'))
