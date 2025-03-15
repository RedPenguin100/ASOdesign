import pandas as pd

from asodesigner.consts import EXPERIMENT_RESULTS


def save_results_organism(df: pd.DataFrame, organism: str, experiment_name: str, algorithm: str):
    csv_path = EXPERIMENT_RESULTS / f'{experiment_name}/{organism}_results/{algorithm}.csv'
    df.to_csv(csv_path, index=False)

def save_results_on_target(df: pd.DataFrame, experiment_name: str, algorithm: str):
    csv_path = EXPERIMENT_RESULTS / f'{experiment_name}/antisense_results/{algorithm}.csv'
    df.to_csv(csv_path, index=False)
