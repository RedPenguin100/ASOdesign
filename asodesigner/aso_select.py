import pandas as pd
import math

from asodesigner.experiment import get_experiments
from asodesigner.read_human_transcriptome import Experiment
from asodesigner.target_finder import get_gfp_second_exp
from asodesigner.util import get_antisense

ANTISENSE_PROPERTIES_CSVS = ['_antisense_fold.csv', 'antisense_melting_temperature.csv',
                             'antisense_nucleotide_properties.csv', 'on_target_fold.csv', 'on_target_energy.csv',
                             '_antisense_self_dimerization_unmodified.csv'
                             ]
OFF_TARGET_PROPERTIES_CSVS = ['gfp_off_targets.csv', 'gfp_hybridization_off_targets.csv']


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


def get_organism_results(organism: str, columns, csv_filename):
    if organism not in ['yeast', 'human']:
        raise ValueError("organism must be yeast or human")

    results_folder = get_results_folder(organism=organism)
    df = pd.read_csv(f"{results_folder}/{experiment.name}{csv_filename}")

    for column in columns:
        df[f"{column}_{organism}"] = df[column]
    df = df.drop(columns, axis=1)
    return df


def load_all_features(experiment: Experiment):
    tables = []
    for csv_path in ANTISENSE_PROPERTIES_CSVS:
        df = pd.read_csv(f"antisense_results/{experiment.name}{csv_path}")
        tables.append(df)

    for csv_filename in OFF_TARGET_PROPERTIES_CSVS:
        if 'gfp_off_targets' in csv_filename:
            columns = ['0_matches', '1_matches', '2_matches', '3_matches']
        elif 'hybridization_off_targets' in csv_filename:
            columns = ['total_hybridization_candidates', 'total_hybridization_energy', 'total_hybridization_max_sum',
                       'total_hybridization_binary_sum']
        else:
            raise ValueError(f"CSV filename {csv_filename} not recognized")

        for organism in ['yeast', 'human']:
            df = get_organism_results(organism=organism, columns=columns, csv_filename=csv_filename)
            tables.append(df)

    merged_df = tables[0]
    for table in tables[1:]:
        merged_df = pd.merge(merged_df, table, on=['sense_start', 'sense_length'])

    print('All columns: ', merged_df.columns)
    all_asos = []
    region_names = []
    for row in merged_df.itertuples():
        i, l = row.sense_start, row.sense_length
        all_asos.append(get_antisense(experiment.target_sequence[i:i + l]))

        if row.sense_start < 100:
            region_names.append('start')
        elif (row.sense_start + row.sense_length) > (
                len(experiment.target_sequence) - 150) and (
                (row.sense_start + row.sense_length) < len(experiment.target_sequence) - 50):
            region_names.append('end_inside')
        elif (row.sense_start + row.sense_length) > (len(experiment.target_sequence) - 50):
            region_names.append('end_outside')
        else:
            region_names.append('middle')

    merged_df['antisense'] = all_asos
    merged_df['region_name'] = region_names

    all_asos_unique = set(all_asos)

    duplicates = len(all_asos) - len(all_asos_unique)
    print(f"Eliminated {duplicates} dups")
    print(all_asos_unique)

    return merged_df


if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)

    experiment_names = ['Second', 'SecondScrambled']
    experiments = get_experiments(experiment_names)

    for experiment in experiments:
        merged_df = load_all_features(experiment)
        # Per Tamir notes
        merged_df = merged_df[merged_df['0_matches_human'] == 0]
        merged_df = merged_df[merged_df['0_matches_yeast'] == 0]

        if experiment.name == 'Second':
            most_important_columns = ['sense_start', 'sense_length', 'on_target_energy_max',
                                      'on_target_fold_openness_normalized', 'total_hybridization_max_sum_human',
                                      'total_hybridization_max_sum_yeast', 'mfe']

            merged_df = merged_df.sort_values(
                by=['on_target_energy_max', 'on_target_fold_openness_normalized', 'total_hybridization_max_sum_human',
                    'total_hybridization_max_sum_yeast', 'mfe'], ascending=[True, False, True, True, False],
                inplace=False)

            irrelevant_columns = ['aso_aso_delta_g', 'aso_aso_delta_h', 'aso_aso_delta_s', 'aso_aso_structure_found',
                                  '2_matches_human', '3_matches_human', '2_matches_yeast', '3_matches_yeast']
            merged_df = merged_df.drop(irrelevant_columns, axis=1)
            reordered_columns = most_important_columns + [col for col in merged_df.columns if
                                                          col not in most_important_columns]
            merged_df = merged_df[reordered_columns]
        elif experiment.name == 'SecondScrambled':
            most_important_columns = ['sense_start', 'sense_length', 'on_target_energy_max',
                                      'on_target_fold_openness_normalized', 'total_hybridization_max_sum_human',
                                      'total_hybridization_max_sum_yeast', 'mfe']

            # On scrambled we would like to miss the target, at least not aim for it
            merged_df = merged_df.sort_values(
                by=['on_target_energy_max', 'on_target_fold_openness_normalized', 'total_hybridization_max_sum_human',
                    'total_hybridization_max_sum_yeast', 'mfe'], ascending=[False, True, True, True, False],
                inplace=False)

            irrelevant_columns = ['aso_aso_delta_g', 'aso_aso_delta_h', 'aso_aso_delta_s', 'aso_aso_structure_found',
                                  '2_matches_human', '3_matches_human', '2_matches_yeast', '3_matches_yeast']
            merged_df = merged_df.drop(irrelevant_columns, axis=1)
            reordered_columns = most_important_columns + [col for col in merged_df.columns if
                                                          col not in most_important_columns]
            merged_df = merged_df[reordered_columns]

        merged_df.to_csv(f'{experiment.name}_all_features.csv', index=False)
