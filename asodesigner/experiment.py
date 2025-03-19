import math
from typing import List, Dict

from asodesigner.consts import EXPERIMENT_RESULTS
from asodesigner.target_finder import get_gfp_second_exp, get_degron_and_gap_third_exp, \
    get_degron_gfp_scrambled_third_exp, get_3utr_gfp
from asodesigner.util import get_antisense


class Experiment:
    def __init__(self):
        self.name = None
        self.target_sequence = None
        self.l_values = None
        self.aso_template = None
        self.gc_content_filter = (-math.inf, math.inf)

    def get_aso_template(self):
        if self.aso_template is None:
            return self.target_sequence
        return self.aso_template


DEFAULT_LENGTHS_UNMODIFIED = [16, 17, 18, 19, 20, 21, 22]


# do not use outside this file
def _get_experiments_dict() -> Dict[str, Experiment]:
    name_to_experiment = {}

    second_scrambled = Experiment()
    second_scrambled.target_sequence = get_gfp_second_exp()
    second_scrambled.aso_template = get_antisense(get_gfp_second_exp())
    second_scrambled.name = 'SecondScrambled'
    second_scrambled.l_values = DEFAULT_LENGTHS_UNMODIFIED
    name_to_experiment[second_scrambled.name] = second_scrambled

    second = Experiment()
    second.target_sequence = get_gfp_second_exp()
    second.name = 'Second'
    second.l_values = DEFAULT_LENGTHS_UNMODIFIED
    name_to_experiment[second.name] = second

    third_degron = Experiment()
    third_degron.target_sequence = get_degron_and_gap_third_exp()
    third_degron.name = 'ThirdDegron'
    third_degron.l_values = DEFAULT_LENGTHS_UNMODIFIED
    name_to_experiment[third_degron.name] = third_degron

    third_gfp_degron_scrambled = Experiment()
    third_gfp_degron_scrambled.target_sequence = get_degron_gfp_scrambled_third_exp()
    # TODO: come with a better idea than aso_template
    third_gfp_degron_scrambled.aso_template = get_antisense(get_degron_gfp_scrambled_third_exp())
    third_gfp_degron_scrambled.name = 'ThirdScrambled'
    third_gfp_degron_scrambled.l_values = DEFAULT_LENGTHS_UNMODIFIED
    name_to_experiment[third_gfp_degron_scrambled.name] = third_gfp_degron_scrambled

    # On target 3UTR
    fourth = Experiment()
    fourth.target_sequence = get_3utr_gfp()
    fourth.name = 'Fourth'
    fourth.l_values = DEFAULT_LENGTHS_UNMODIFIED
    name_to_experiment[fourth.name] = fourth

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
