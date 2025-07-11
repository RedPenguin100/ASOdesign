import pandas as pd
import os

def load_all_features(filenames=None):
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


def save_my_features(df, feature_name):
    # this function receives the df of your feature(s) (with indices related to the rows please) and saves it to our feature directory
    # to use this add - from scripts.features.feature_extraction import save_my_features
    project_root = find_project_root('scripts')
    feature_dir = os.path.join(project_root, 'scripts', 'features')
    os.makedirs(feature_dir, exist_ok=True)
    file_path = os.path.join(feature_dir, f'{feature_name}.csv')
    df.to_csv(file_path, index=False)


if __name__ == "__main__":
    df1 = pd.DataFrame({'index': [1, 2, 3], 'A': ['a1', 'a2', 'a3']})
    df2 = pd.DataFrame({'index': [2, 3, 4], 'B': ['b2', 'b3', 'b4']})
    df3 = pd.DataFrame({'index': [1, 2, 4], 'C': ['c1', 'c2', 'c4']})

    dfs = [df1, df2, df3]

    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='index', how='outer')

    print(merged_df)