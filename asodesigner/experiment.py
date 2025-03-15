import os
from typing import List, Dict

from asodesigner.consts import EXPERIMENT_RESULTS
from asodesigner.target_finder import get_gfp_second_exp
from asodesigner.util import get_antisense


class Experiment:
    def __init__(self):
        self.name = None
        self.target_sequence = None
        self.l_values = None


# do not use outside this file
def _get_experiments_dict() -> Dict[str, Experiment]:
    name_to_experiment = {}

    second_scrambled = Experiment()
    second_scrambled.target_sequence = get_antisense(get_gfp_second_exp())
    second_scrambled.name = 'SecondScrambled'
    second_scrambled.l_values = [16, 17, 18, 19, 20, 21, 22]
    name_to_experiment[second_scrambled.name] = second_scrambled

    second = Experiment()
    second.target_sequence = get_gfp_second_exp()
    second.name = 'Second'
    second.l_values = [16, 17, 18, 19, 20, 21, 22]
    name_to_experiment[second.name] = second
    return name_to_experiment


def maybe_create_experiment_folders(experiment_name: str):
    """
    maybe in function name because we don't want to fail if experiment exists
    :param experiment_name:
    """
    experiment_path = EXPERIMENT_RESULTS / experiment_name
    yeast_results = experiment_path / 'yeast_results'
    human_results = experiment_path / 'human_results'
    antisense_results = experiment_path / 'antisense_results'

    experiment_path.mkdir(exist_ok=True)
    yeast_results.mkdir(exist_ok=True)
    human_results.mkdir(exist_ok=True)
    antisense_results.mkdir(exist_ok=True)

    print(f"Experiment {experiment_name} created successfully.")


def get_experiment(name: str) -> Experiment:
    name_to_experiment = _get_experiments_dict()
    return name_to_experiment[name]


def get_experiments(names) -> List[Experiment]:
    name_to_experiment = _get_experiments_dict()

    if names is None:
        return list(name_to_experiment.values())

    experiments = []
    for name in names:
        experiments.append(name_to_experiment[name])

    return experiments
