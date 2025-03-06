import pandas as pd
import math

from Bio import SeqIO

from asodesigner.target_finder import get_gfp_first_exp
from asodesigner.util import get_antisense

ANTISENSE_PROPERTIES_CSVS = ['antisense_fold.csv', 'antisense_melting_temperature.csv',
                             'antisense_nucleotide_properties.csv',
                             ]
ENVIRONMENT_PROPERTIES_CSVS = ['gfp_off_targets.csv']


def filter_gc_content(df, min_content=-math.inf, max_content=math.inf):
    return df[(min_content <= df['gc_content']) & (df['gc_content'] <= max_content)]


def filter_weakly_folded(df, min_mfe=0.):
    return df[df['mfe'] >= min_mfe]


def filter_off_targets(df, max0=0, max1=math.inf, max2=math.inf, max3=math.inf):
    return df[
        (df['0_matches'] <= max0) & (df['1_matches'] <= max1) & (df['2_matches'] <= max2) & (df['3_matches'] <= max3)]


def filter_length(df, min_length=15, max_length=22):
    return df[(df['sense_length'] >= min_length) & (df['sense_length'] <= max_length)]


def filter_melting_temperature(df, min_temp=-math.inf, max_temp=math.inf):
    return df[(df['melting_temperature'] >= min_temp) & (df['melting_temperature'] <= max_temp)]


def get_results_folder(organism: str) -> str:
    if organism == "yeast":
        return "yeast_results"
    elif organism == "human":
        return "human_results"
    else:
        raise ValueError("organism must be yeast or human")


if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)

    target_seq = get_gfp_first_exp()
    exp_name = 'First'

    tables = []
    for csv_path in ANTISENSE_PROPERTIES_CSVS:
        df = pd.read_csv(f"antisense_results/{exp_name}{csv_path}")
        tables.append(df)

    for csv_filename in ENVIRONMENT_PROPERTIES_CSVS:
        results_folder = get_results_folder(organism='yeast')
        df = pd.read_csv(f"{results_folder}/{exp_name}{csv_filename}")
        df['0_matches_yeast'] = df['0_matches']
        df['1_matches_yeast'] = df['1_matches']
        df['2_matches_yeast'] = df['2_matches']
        df['3_matches_yeast'] = df['3_matches']
        df = df.drop(['0_matches', '1_matches', '2_matches', '3_matches'], axis=1)
        tables.append(df)

        results_folder = get_results_folder(organism='human')
        df = pd.read_csv(f"{results_folder}/{exp_name}{csv_filename}")
        df['0_matches_human'] = df['0_matches']
        df['1_matches_human'] = df['1_matches']
        df['2_matches_human'] = df['2_matches']
        df['3_matches_human'] = df['3_matches']
        df = df.drop(['0_matches', '1_matches', '2_matches', '3_matches'], axis=1)

        tables.append(df)

    merged_df = tables[0]
    for table in tables[1:]:
        merged_df = pd.merge(merged_df, table, on=['sense_start', 'sense_length'])

    print('All columns: ', merged_df.columns)

    merged_df = merged_df.drop(
        ['structure', '3_matches_yeast', '3_matches_human', '2_matches_yeast', '2_matches_human'], axis=1)

    merged_df = merged_df.sort_values(by=['0_matches_human', '0_matches_yeast', 'mfe', 'sense_start'], ascending=[True, True, False, True])
    merged_df = merged_df[(merged_df['0_matches_human'] == 0) & (merged_df['0_matches_yeast'] == 0)]


    # merged_df = merged_df[(merged_df['sense_start']  > (len(target_seq) - 200)) & (merged_df['sense_start'] < (len(target_seq) - 100))]
    # merged_df = merged_df[(merged_df['sense_start'] < 250) & (merged_df['sense_start'] > 100)]
    merged_df = merged_df[(merged_df['sense_start'] > 790)]


    merged_df = filter_gc_content(merged_df, 0.45, 0.55)
    merged_df = filter_length(merged_df, 20, 20)
    # merged_df = filter_weakly_folded(merged_df)

    # merged_df = filter_weakly_folded(merged_df)
    # merged_df = filter_gc_content(merged_df, 0.45, 0.55)
    # merged_df = filter_melting_temperature(merged_df, 45, 55)
    # merged_df = filter_off_targets(merged_df, max0=0, max1=0, max2=0)
    # df_14 = filter_length(merged_df, length=14)
    # df_15 = filter_length(merged_df, length=15)
    # df_16 = filter_length(merged_df, length=16)
    # print(df_14)
    # print(df_15)
    # print(df_16)

    all_asos = []
    for row in merged_df.itertuples():
        i, l = row.sense_start, row.sense_length
        all_asos.append(get_antisense(target_seq[i:i + l]))
    merged_df['antisense'] = all_asos
    print(merged_df.head(30))

    all_asos_unique = set(all_asos)

    duplicates = len(all_asos) - len(all_asos_unique)

    print(f"Eliminated {duplicates} dups")
    print(all_asos_unique)
