import ViennaRNA as RNA
import numpy as np

import os
import re
import random
import subprocess
from Bio.Seq import Seq
from typing import Dict, List

from asodesigner.consts import RISEARCH1_BINARY_PATH


def get_weighted_energy(target_start, l, step_size, energies, window_size):
    target_end = target_start + l
    proxy = np.zeros(l, dtype=np.float64)
    window_residual = window_size % step_size

    for j in range(target_start, target_end):
        proxy_index = j - target_start
        j_residual = j % step_size
        j_energy = j // step_size

        if j_energy < (window_size // step_size):
            proxy[proxy_index] = np.mean(energies[0:(j_energy + 1)], dtype=np.float64)
        else:
            if j_residual <= (window_residual - 1):
                proxy[proxy_index] = np.mean(
                    [energies[j_energy - 2], energies[j_energy - 1], energies[j_energy]], dtype=np.float64)
            else:
                proxy[proxy_index] = np.mean(
                    [energies[j_energy - 1], energies[j_energy]], dtype=np.float64)
    return np.mean(proxy, dtype=np.float64)


def calculate_energies(target_seq, step_size, window_size):
    energies = np.empty(len(target_seq) // step_size + 1)
    for i in range(0, len(target_seq) - window_size + 1, step_size):
        _, mfe = RNA.fold(target_seq[i: i + window_size])
        energies[i // step_size] = mfe
    return energies


'''
note: I used hashing here to be able to run multiple sequences (triggers in my case) in parallel (for multiple triggers) without conflicts with the files,
you can remove it if it's not an issue for you 
note: in my case I first preformed reverse_complement_rna to the trigger because my goal was to find sequences that might undesirably bind to the
trigger binding site of the toehold, adjust it for your particular case
note: play with the d and s parameters I passed here, and with other parameters that might be relevant for your usecase
'''


# in my case I was only interested in the energy scores and didn't care about the actual sequences, this is the parsing I used in case it helps you
def get_mfe_scores(result: str) -> List[List[float]]:
    mfe_results = []

    for gene_result in result.split("\n\nquery trigger")[1:]:
        stripped_result = gene_result.strip()
        regex_results = re.findall("Free energy \[kcal/mol\]: [0-9-.]+ ", stripped_result)
        mfe_results.append(
            [float(regex_result.replace('Free energy [kcal/mol]: ', '').strip()) for regex_result in regex_results])

    return mfe_results


def get_trigger_mfe_scores_by_risearch(trigger: str, name_to_sequence: Dict[str, str]) -> str:
    # logger.info(f"started evaluating a single trigger: {trigger}")
    hash = random.getrandbits(64)

    with open(f"target-{hash}.fa", "w") as f:
        for name, sequence in name_to_sequence.items():
            f.write(">" + str(name) + "\n" + sequence + "\n")

    with open(f"query-{hash}.fa", "w") as f:
        f.write(">trigger" + "\n" + str(Seq(trigger).reverse_complement_rna()) + "\n")

    result = subprocess.check_output(
        [RISEARCH1_BINARY_PATH, "-q", f"query-{hash}.fa", "-t", f"target-{hash}.fa", "-s", "1200", "-d", "30"],
        universal_newlines=True,
        text=True
    )
    os.remove(f"target-{hash}.fa")
    os.remove(f"query-{hash}.fa")
    # logger.info(f"finished evaluating a single trigger: {trigger}")
    return result
