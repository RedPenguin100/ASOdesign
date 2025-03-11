from typing import List

from asodesigner.target_finder import get_gfp_second_exp
from asodesigner.util import get_antisense


class Experiment:
    def __init__(self):
        self.name = None
        self.target_sequence = None
        self.l_values = None


def get_experiments(names) -> List[Experiment]:
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

    if names is None:
        return list(name_to_experiment.values())

    experiments = []
    for name in names:
        experiments.append(name_to_experiment[name])

    return experiments
