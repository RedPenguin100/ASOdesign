import math
import bisect

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
from more_itertools import divide
from tqdm import tqdm
import json

from asodesigner.consts import YEAST_FASTA_PATH, YEAST_FIVE_PRIME_UTR, YEAST_THREE_PRIME_UTR, GFP1_PATH, \
    YEAST_GFF_DB_PATH
from asodesigner.consts import YEAST_GFF_PATH
from fuzzysearch import find_near_matches

import pandas as pd

from asodesigner.read_utils import process_hybridization
from asodesigner.fold import get_weighted_energy, calculate_energies, get_trigger_mfe_scores_by_risearch, get_mfe_scores
from asodesigner.target_finder import get_gfp_first_exp, get_gfp_second_exp
from asodesigner.timer import Timer
from asodesigner.util import get_longer_string


def cond_print(text, verbose=False):
    if verbose:
        print(text)


class LocusInfo:
    def __init__(self):
        self.exons = []
        self.introns = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.exon_concat = None


def get_locus_to_data_dict_alternative():
    db_path = Path(YEAST_GFF_DB_PATH)

    if not db_path.exists():
        with Timer() as t:
            db = gffutils.create_db(str(YEAST_GFF_PATH), dbfn=str(db_path), force=True, keep_order=True,
                                    merge_strategy='merge', sort_attribute_values=True)
            # db.update(list(db.create_introns())) # Uncomment to create introns
        print(f"DB create took: {t.elapsed_time}s")
    else:
        print("Opening DB")
        db = gffutils.FeatureDB(str(db_path))

    fasta_dict = SeqIO.to_dict(SeqIO.parse(str(YEAST_FASTA_PATH), 'fasta'))
    locus_to_data = dict()
    locus_to_strand = dict()

    for feature in db.features_of_type(('CDS', 'intron')):
        chrom = feature.seqid
        if chrom == 'NC_001224.1':  # skip mitochondrial DNA
            continue

        locus_tags = feature.attributes['locus_tag']
        if len(locus_tags) != 1:
            raise ValueError(f"Multiple locuses: {locus_tags}")
        locus_tag = locus_tags[0]

        seq = fasta_dict[chrom].seq[feature.start - 1: feature.end]

        if feature.strand == '-':
            seq = seq.reverse_complement()
        seq = seq.upper()

        if feature.featuretype == 'CDS':
            cds = feature

            if locus_tag not in locus_to_data:
                locus_info = LocusInfo()
                locus_info.exons = [(cds.start, seq)]
                locus_to_data[locus_tag] = locus_info

                locus_to_strand[locus_tag] = cds.strand
            else:
                bisect.insort(locus_to_data[locus_tag].exons, (cds.start, seq))
        elif feature.featuretype == 'intron':
            intron = feature

            if locus_tag not in locus_to_data:
                locus_info = LocusInfo()
                locus_info.introns = [(feature.start, seq)]
                locus_to_data[locus_tag] = locus_info

                locus_to_strand[locus_tag] = intron.strand
            else:
                bisect.insort(locus_to_data[locus_tag].introns, (intron.start, seq))

    for locus_tag in locus_to_data:
        locus_info = locus_to_data[locus_tag]
        if locus_to_strand[locus_tag] == '-':
            locus_info.exons.reverse()
            locus_info.introns.reverse()

        locus_info = locus_to_data[locus_tag]
        locus_info.exons = [element for _, element in locus_info.exons]
        locus_info.introns = [element for _, element in locus_info.introns]
        locus_info.exon_concat = "".join(str(seq) for seq in locus_info.exons)

    return locus_to_data


def get_locus_to_data_dict():
    locus_to_data = dict()

    with Timer() as t:
        with open(YEAST_FASTA_PATH, 'r') as fasta_handle, open(YEAST_GFF_PATH, 'r') as gff_handle:
            for record in GFF.parse(gff_handle, base_dict=SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))):
                if "NC_001224.1" == record.id:  # Mitochondrial DNA, skip
                    continue

                cond_print(f"> {record.id}")

                cond_print(record.features)
                for feature in record.features:
                    cond_print(f"> {feature.type}")
                    if feature.type == 'gene':
                        cond_print(f" > {feature}")
                        cond_print(f" > {feature.id}")
                        for sub_feature in feature.sub_features:
                            for sub_sub_feature in sub_feature.sub_features:
                                cond_print(f">> {sub_sub_feature}")
                                if sub_sub_feature.type == 'CDS':
                                    cond_print(f"  CDS ID: {sub_sub_feature.id}")
                                    cond_print(f"  Location: {sub_sub_feature.location}")
                                    cond_print(f"  Parent Gene: {sub_sub_feature.qualifiers.get('Parent')}")
                                    if len(sub_sub_feature.qualifiers['locus_tag']) != 1:
                                        raise ValueError(f"More than 1 locus tag on CDS ID: {sub_sub_feature.id}")
                                    locus_tag = sub_sub_feature.qualifiers['locus_tag'][0]
                                    seq = sub_sub_feature.extract(record.seq)

                                    locus_to_data.setdefault(locus_tag, LocusInfo()).exons.append(seq.upper())

    print(f"Took: {t.elapsed_time}s")

    return locus_to_data


def process_sense_one_locus(args):
    lvalues, gfp_seq, locus_data_chunk = args
    matches_per_distance = [0, 0, 0, 0]

    df = pd.DataFrame(
        columns=['sense_start', 'sense_length', '0_matches', '1_matches', '2_matches', '3_matches']).set_index(
        ['sense_start', 'sense_length'])

    count = 0
    split_size = 100

    for locus_data in locus_data_chunk:
        locus_name, locus_info = locus_data

        for l in lvalues:
            for i in range(len(gfp_seq) - l + 1):
                sense = gfp_seq[i:i + l]

                matches = find_near_matches(sense, locus_info.full_mrna, max_insertions=0, max_deletions=0, max_l_dist=3)
                for match in matches:
                    matches_per_distance[match.dist] += 1

                if (i, l) not in df.index:
                    df.loc[len(df)] = [i, l, matches_per_distance[0], matches_per_distance[1], matches_per_distance[2],
                                       matches_per_distance[3]]

        count += 1
        if count == split_size:
            df = df.groupby(['sense_start', 'sense_length'], as_index=False).sum()
            count = 0

    return df


def process_sense_locus_off_target(args):
    i, l, target_seq, locus_to_data = args
    sense = target_seq[i:i + l]
    matches_per_distance = [0, 0, 0, 0]

    for locus_tag, locus_info in locus_to_data.items():
        matches = find_near_matches(sense, locus_info, max_insertions=0, max_deletions=0, max_l_dist=3)
        for match in matches:
            matches_per_distance[match.dist] += 1

    # Return a tuple containing the starting index, current l, and match counts
    return (i, l, matches_per_distance[0],
            matches_per_distance[1], matches_per_distance[2], matches_per_distance[3])


def process_fold_single_mrna(args):
    locus_tag, locus_info, step_size, window_size = args
    energies = calculate_energies(locus_info.full_mrna, step_size=step_size, window_size=window_size)
    return locus_tag, energies


def process_locus_fold_off_target(locus_to_data):
    window_size = 40
    step_size = 15

    results = dict()
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_fold_single_mrna, (arg[0], arg[1], step_size, window_size)) for arg in
                   locus_to_data.items()]

        for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
            locus, energies = future.result()
            results[locus] = energies.tolist()

    with open(f'yeast_results/Firstgfp_fold_off_targets_window_{window_size}_step_{step_size}.csv', 'w') as f:
        json.dump(results, f, indent=4)

def load_three_prime_utr(locus_to_data):
    for record in SeqIO.parse(YEAST_THREE_PRIME_UTR, "fasta"):
        locus = record.name.split('_')[4]
        if locus in locus_to_data:
            locus_info = locus_to_data[locus]

            locus_info.three_prime_utr = get_longer_string(record.seq.upper()[:-1],
                                                           locus_info.three_prime_utr)


def load_five_prime_utr(locus_to_data):
    for record in SeqIO.parse(YEAST_FIVE_PRIME_UTR, "fasta"):
        locus = record.name.split('_')[4]
        if locus in locus_to_data:
            locus_info = locus_to_data[locus]

            locus_info.five_prime_utr = get_longer_string(record.seq.upper()[:-1],
                                                          locus_info.five_prime_utr)


if __name__ == "__main__":
    locus_to_data = get_locus_to_data_dict_alternative()

    load_five_prime_utr(locus_to_data)
    load_three_prime_utr(locus_to_data)
    for locus_name, locus_info in locus_to_data.items():
        locus_info.full_mrna = f"{locus_info.five_prime_utr}{locus_info.exon_concat}{locus_info.three_prime_utr}"

    target_seq = get_gfp_second_exp()
    experiment_name = 'Second'

    tasks = []
    l_values = [15, 16, 17, 18, 19, 20, 21, 22]

    max_iterations = math.inf

    simple_locus_to_data = dict()
    for i in range(len(locus_to_data)):
        locus_name, locus_info = list(locus_to_data.items())[i]
        simple_locus_to_data[locus_name] = locus_info.full_mrna

    for l in l_values:
        for i in range(min(len(target_seq) - l + 1, max_iterations)):
            tasks.append((i, l, target_seq, simple_locus_to_data))

    # process_locus_fold_off_target(locus_to_data)

    # with Timer() as t:
    #     results = []
    #     with ProcessPoolExecutor() as executor:
    #         futures = [executor.submit(process_hybridization, arg) for arg in tasks]
    #
    #         for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
    #             results.append(future.result())
    #
    # columns = ['sense_start', 'sense_length', 'total_hybridization_candidates', 'total_off_target_energy',
    #            'total_max_sum', 'total_binary_sum']
    # df = pd.DataFrame(results, columns=columns)
    # df.to_csv(f'yeast_results/{experiment_name}hybridization_candidates3.csv', index=False)

    columns = ['sense_start', 'sense_length', '0_matches', '1_matches', '2_matches', '3_matches']

    results = []
    with Timer() as t:
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_sense_locus_off_target, arg) for arg in tasks]
            for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):

                results.append(future.result())

    df = pd.DataFrame(results, columns=columns)

    aggregated_df = df.groupby(['sense_start', 'sense_length'], as_index=False).sum()


    df = aggregated_df
    df.to_csv('yeast_results/gfp_off_targets.csv', index=False)


    print(f"Time took to find: {t.elapsed_time}s")

    print(df)
