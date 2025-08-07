import math
import random
from typing import List, Dict

import pandas as pd

from asodesigner.consts import EXPERIMENT_RESULTS, DATA_PATH
from asodesigner.random_util import generate_random_dna
from asodesigner.target_finder import get_gfp_second_exp, get_degron_and_gap_third_exp, \
    get_degron_gfp_scrambled_third_exp, get_3utr_gfp, get_extended_gfp, iterate_template, iterate_template_antisense, \
    get_angtpl2, get_bcl2_patent
from asodesigner.util import get_antisense


class ExperimentSetting:
    TEMPLATE = "TEMPLATE"
    GENERATED = "GENERATED"

class Experiment:
    def __init__(self):
        self.name = None
        self.target_sequence = None
        self.l_values = None
        self.aso_template = None
        self.gc_content_filter = (-math.inf, math.inf)
        self.setting = ExperimentSetting.TEMPLATE  # can either be TEMPLATE or GENERATED
        self.generated_list = None

    def get_aso_template(self):
        if self.aso_template is None:
            return self.target_sequence
        return self.aso_template

    def get_aso_sense_iterator(self):
        if self.setting == ExperimentSetting.TEMPLATE:
            # NOTE: only yield from will defer the iterator properly
            yield from iterate_template(self.get_aso_template(), self.l_values)
        elif self.setting == ExperimentSetting.GENERATED:
            for i in range(len(self.generated_list)):
                sense = self.generated_list[i]
                yield i, len(sense), sense
        else:
            raise ValueError(f"Unknown experiment setting: {self.setting}")

    def get_aso_antisense_iterator(self):
        if self.setting == ExperimentSetting.TEMPLATE:
            # NOTE: only yield from will defer the iterator properly
            yield from iterate_template_antisense(self.get_aso_template(), self.l_values)
        elif self.setting == ExperimentSetting.GENERATED:
            for i in range(len(self.generated_list)):
                antisense = get_antisense(self.generated_list[i])
                yield i, len(antisense), antisense
        else:
            raise ValueError(f"Unknown experiment setting: {self.setting}")

    def get_aso_sense_by_index(self, idx, length):
        if self.setting == ExperimentSetting.TEMPLATE:
            return self.get_aso_template()[idx:idx + length]
        elif self.setting == ExperimentSetting.GENERATED:
            return self.generated_list[idx]
        raise ValueError(f"Unknown experiment setting: {self.setting}")

    def get_aso_antisense_by_index(self, idx, length):
        if self.setting == ExperimentSetting.TEMPLATE:
            return get_antisense(self.get_aso_template()[idx:idx + length])
        elif self.setting == ExperimentSetting.GENERATED:
            return get_antisense(self.generated_list[idx])
        raise ValueError(f"Unknown experiment setting: {self.setting}")


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

    # Finally
    entire = Experiment()
    entire.target_sequence = get_extended_gfp()
    entire.name = 'Entire'
    entire.l_values = DEFAULT_LENGTHS_UNMODIFIED
    name_to_experiment[entire.name] = entire

    # NOTE: I accidentally generated the populated sense values with antisense values.
    # everything is ok for the experiment, but the fold property does not necessarily hold, and
    # needs to be verified when picking values.
    # future versions need to use `[get_antisense(antisense) for antisense in generate_random_dna(...)]`
    entire_scrambled = Experiment()
    entire_scrambled.target_sequence = get_extended_gfp()
    entire_scrambled.name = 'EntireScrambled'
    entire_scrambled.setting = ExperimentSetting.GENERATED
    entire_scrambled.l_values = DEFAULT_LENGTHS_UNMODIFIED

    random.seed(0)
    entire_scrambled.generated_list = []
    for l in entire_scrambled.l_values:
        entire_scrambled.generated_list.extend(generate_random_dna(l, attempts=400))
    name_to_experiment[entire_scrambled.name] = entire_scrambled

    angptl2_exp = Experiment()
    angptl2_exp.target_sequence = get_angtpl2()
    angptl2_exp.name = 'Angptl2'
    angptl2_exp.setting = ExperimentSetting.GENERATED
    angptl2_exp.l_values = DEFAULT_LENGTHS_UNMODIFIED  # TODO: remove this
    angptl2_exp.generated_list = [get_antisense(sense) for sense in
                                  list(pd.read_csv(DATA_PATH / 'ANGPTL2_antisense.csv')['Sequence'])]
    name_to_experiment[angptl2_exp.name] = angptl2_exp

    positive_control = Experiment()
    positive_control.target_sequence = get_extended_gfp()
    positive_control.name = 'EntirePositiveControl'
    positive_control.setting = ExperimentSetting.GENERATED
    positive_control.l_values = DEFAULT_LENGTHS_UNMODIFIED  # TODO: remove this

    positive_control.generated_list = [
        # Fix for ASO5 in this article https://ars.els-cdn.com/content/image/1-s2.0-S2162253120300536-mmc1.pdf
        get_antisense('TTGCCGGTGGTGCAGATGAA'),
        # Original ASO5
        get_antisense('TTGCCGGTGGTGCAGATAAA'),

        get_antisense('TGTGGCGGATCTTGAAGTTC'),
        get_antisense('CTGCTGGTAGTGGTCGGCGA'),
        get_antisense('GCGGACTGGGTGCTCAGGTA'),
        get_antisense('ACGATGGTCCTTCTTGTGAC'),

    ]

    name_to_experiment[positive_control.name] = positive_control


    bcl2_patent_exp = Experiment()
    bcl2_patent_exp.target_sequence = get_bcl2_patent()
    bcl2_patent_exp.name = 'BCL2Patent'
    bcl2_patent_exp.setting = ExperimentSetting.GENERATED
    bcl2_patent_exp.l_values = DEFAULT_LENGTHS_UNMODIFIED # TODO: remove this
    bcl2_patent_exp.generated_list = [get_antisense(sense) for sense in list(pd.read_csv(DATA_PATH / 'experiments_bcl2.csv')['Sequence'])]
    name_to_experiment[bcl2_patent_exp.name] = bcl2_patent_exp



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
