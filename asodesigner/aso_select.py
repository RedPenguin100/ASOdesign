import pandas as pd
import math

from asodesigner.consts import EXPERIMENT_RESULTS
from asodesigner.experiment import get_experiments, Experiment
from asodesigner.features.feature_names import SENSE_LENGTH, SENSE_START
from asodesigner.util import get_antisense


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


def get_explicit_antisense(merged_df, experiment):
    all_asos = []
    for row in merged_df.itertuples():
        i, l = row.sense_start, row.sense_length
        all_asos.append(experiment.get_aso_antisense_by_index(idx=i, length=l))
    return all_asos


def mark_regions(merged_df, experiment):
    region_names = []
    aso_template = experiment.get_aso_template()

    if experiment.name == "Second":
        for row in merged_df.itertuples():
            if row.sense_start < 100:
                region_names.append('start')
            elif (row.sense_start + row.sense_length) > (
                    len(aso_template) - 150) and (
                    (row.sense_start + row.sense_length) < len(aso_template) - 50):
                region_names.append('end_inside')
            elif (row.sense_start + row.sense_length) > (len(aso_template) - 50):
                region_names.append('end_outside')
            else:
                region_names.append('middle')
    elif experiment.name == 'Entire':
        for row in merged_df.itertuples():
            if row.sense_start < 100:
                region_names.append('GFP_start')
            elif (row.sense_start + row.sense_length) < 500:
                region_names.append('GFP_middle')
            elif (row.sense_start + row.sense_length) < 717:
                region_names.append('GFP_end')
            elif (row.sense_start + row.sense_length) < 717 + 126:
                region_names.append('Degron')
            elif (row.sense_start + row.sense_length) > (len(aso_template) - 132):
                region_names.append('Terminator')
            elif (row.sense_start + row.sense_length) > (len(aso_template) - 132 - 77):
                region_names.append('Gap-(LTR/Term)')
            elif (row.sense_start + row.sense_length) > (len(aso_template) - 132 - 77 - 234):
                region_names.append('3LTR')
            else:
                region_names.append('Other')
    else:
        for _ in merged_df.itertuples():
            region_names.append('NA')

    return region_names


def load_all_features(experiment: Experiment):
    current_experiment_results = EXPERIMENT_RESULTS / experiment.name
    antisense_results = current_experiment_results / 'antisense_results'

    tables = []
    for csv_path in antisense_results.glob('*.csv'):
        df = pd.read_csv(csv_path)
        tables.append(df)

    yeast_results = current_experiment_results / 'yeast_results'
    human_results = current_experiment_results / 'human_results'

    KEY_COLUMNS = [SENSE_START, SENSE_LENGTH]

    for csv_path in yeast_results.glob('*.csv'):
        df = pd.read_csv(csv_path)
        rename_map = {col: f'{col}_yeast' for col in df.columns if col not in KEY_COLUMNS}
        df.rename(columns=rename_map, inplace=True)
        tables.append(df)

    for csv_path in human_results.glob('*.csv'):
        df = pd.read_csv(csv_path)
        rename_map = {col: f'{col}_human' for col in df.columns if col not in KEY_COLUMNS}
        df.rename(columns=rename_map, inplace=True)
        tables.append(df)

    merged_df = tables[0]
    for table in tables[1:]:
        merged_df = pd.merge(merged_df, table, on=KEY_COLUMNS)

    print('All columns: ', merged_df.columns)

    all_asos = get_explicit_antisense(merged_df, experiment)
    region_names = mark_regions(merged_df, experiment)

    merged_df['antisense'] = all_asos
    merged_df['region_name'] = region_names

    all_asos_unique = set(all_asos)

    duplicates = len(all_asos) - len(all_asos_unique)
    print(f"Eliminated {duplicates} dups")
    print(all_asos_unique)

    return merged_df


irrelevant_columns = ['aso_aso_delta_g', 'aso_aso_delta_h', 'aso_aso_delta_s', 'aso_aso_structure_found',
                      '2_matches_human', '3_matches_human', '2_matches_yeast', '3_matches_yeast']


def combine_experiments(experiments):
    experiment_dfs = []
    experiment_names = []

    for experiment in experiments:
        merged_df = load_all_features(experiment)
        # Best way to distinguish experiments in the same table
        merged_df[SENSE_START] = experiment.name + merged_df[SENSE_START].astype(str)
        experiment_dfs.append(merged_df)
        experiment_names.append(experiment.name)
        print(merged_df)


if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)

    experiment_names = ['EntirePositiveControl']
    # experiment_names = ['Second', 'ThirdDegron', 'Fourth']
    experiments = get_experiments(experiment_names)
    # combine_experiments(experiments)

    for experiment in experiments:
        merged_df = load_all_features(experiment)
        # Per Tamir notes
        # merged_df = merged_df[merged_df['0_matches_human'] == 0]
        # merged_df = merged_df[merged_df['0_matches_yeast'] == 0]

        if experiment.name == 'Second' or experiment.name == 'Entire':
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

        elif experiment.name == 'SecondScrambled' or experiment.name == 'EntireScrambled':
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

        experiment_path = EXPERIMENT_RESULTS / experiment.name
        merged_df.to_csv(experiment_path / f'all_features.csv', index=False)
