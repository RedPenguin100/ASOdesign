import pandas as pd
import os

from consts import DATA_PATH_NEW


def load_all_features(filenames=None, light=True):
    # this function loads all the features and the current data, merges all the features according to the index
    # and returns a merged df with all the features and the data
    # missing values would be considered as Nan!

    # to use this add - from scripts.features.feature_extraction import load_all_features
    project_root = find_project_root('scripts')
    feature_dir = os.path.join(project_root, 'scripts', 'features')
    if not filenames:
        filenames = [f for f in os.listdir(feature_dir) if f.endswith('.csv') and not f.startswith('.')]
        filenames.sort()

    if not filenames:
        raise FileNotFoundError(f"No CSV files found in {feature_dir}")

    if light:
        filenames.remove('Smiles.csv')
    print(f"Loading features from: {filenames}")

    dfs = [pd.read_csv(os.path.join(feature_dir, f)) for f in filenames]

    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='index', how='outer')

    return merged_df


def find_project_root(target_folder='scripts'):
    """
    Walk up the directory tree to find the project root,
    which is defined as the folder that contains `scripts/`.
    """
    current_path = os.path.abspath(os.getcwd())
    while True:
        if os.path.isdir(os.path.join(current_path, target_folder)):
            return current_path
        parent = os.path.dirname(current_path)
        if parent == current_path:
            raise FileNotFoundError(f"Could not find folder '{target_folder}' in any parent directory.")
        current_path = parent


def save_feature(df, feature_name):
    project_root = find_project_root('scripts')
    feature_dir = os.path.join(project_root, 'scripts', 'features')
    os.makedirs(feature_dir, exist_ok=True)
    sub_df = df[['index', feature_name]]
    file_path = os.path.join(feature_dir, f'{feature_name}.csv')
    sub_df.to_csv(file_path, index=False)


def read_base_df():
    all_data = pd.read_csv(DATA_PATH_NEW / 'data_asoptimizer_updated.11.7.csv',
                           dtype={'index': int, 'ISIS': int, 'Target_gene': str,
                                  'Cell_line': str, 'Density(cells_per_well)': str,
                                  'Transfection': str, 'ASO_volume(nm)': float, 'Treatment_Period(hours)': float,
                                  'Primer_probe_set': str, 'Sequence': str, 'Modification': str, 'Location': str,
                                  'Chemical_Pattern': str, 'Linkage': str, 'Linage_Location' : str, 'Smiles' : str,
                                  'Inhibition(%)' : float, 'seq_length' : int, 'Canonical Gene Name' : str,
                                  'Cell line organism' : str, 'Transcript' : str, 'Location_in_sequence' : int,
                                  'Location_div_by_length' : float, 'true_length_of_seq' : int, 'mod_scan' : int,
                                  'cell_line_uniform' : str
                                  })
    return all_data


if __name__ == "__main__":
    df1 = pd.DataFrame({'index': [1, 2, 3], 'A': ['a1', 'a2', 'a3']})
    df2 = pd.DataFrame({'index': [2, 3, 4], 'B': ['b2', 'b3', 'b4']})
    df3 = pd.DataFrame({'index': [1, 2, 4], 'C': ['c1', 'c2', 'c4']})

    dfs = [df1, df2, df3]

    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='index', how='outer')

    print(merged_df)
