from asodesigner.fold import get_trigger_mfe_scores_by_risearch

import pandas as pd
from io import StringIO
import os


def dna_to_rna_reverse_complement(seq: str) -> str:
    seq = seq.upper()
    translation_table = str.maketrans("ATGC", "TACG")
    # Translate and reverse
    return seq.translate(translation_table)[::-1]


def parse_risearch_output(output_str: str) -> pd.DataFrame:
    columns = ["trigger", "trigger_start", "trigger_end", "target", "target_start", "target_end", "score", "energy"]
    df = pd.read_csv(StringIO(output_str.strip()), sep="\t", header=None, names=columns)
    return df

def aggregate_off_targets(df: pd.DataFrame) -> pd.DataFrame:
    # Aggregate: sum score (if it's hybridization hits) and take minimum (strongest) energy
    grouped = df.groupby("target").agg({
        "score": "sum",
        "energy": "min"
    }).reset_index()
    return grouped