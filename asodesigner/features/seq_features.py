import re

from Bio.Seq import Seq
from numba import njit
from scipy.stats import entropy

from features.suffix_array import longest_prefix
import numpy as np

import codonbias as cb

from util import get_antisense


@njit
def _is_palindrome(seq: str) -> bool:
    final = ''
    test = seq[::-1]
    for i in range(len(test)):
        if test[i] == 'A':
            final += 'T'
        elif test[i] == 'C':
            final += 'G'
        elif test[i] == 'G':
            final += 'C'
        elif test[i] == 'T':
            final += 'A'
    return final == seq


@njit
def palindromic_fraction(seq: str, l: int) -> float:
    count = 0
    if ~l % 2:
        for n in range(len(seq) - l + 1):
            curr_seq = seq[n:n + l]
            count += _is_palindrome(curr_seq)
    return count / len(seq)


@njit
def homooligo_count(seq: str) -> float:
    seq += '$'
    tot_count = 0
    curr_seq = ''
    n = 0
    while n in range(len(seq) - 1):
        while seq[n] == seq[n + 1]:
            curr_seq += seq[n]
            n += 1
        if len(curr_seq) > 1:
            tot_count += len(curr_seq) + 1
        n += 1
        curr_seq = ''
    return tot_count / len(seq)


def compute_ENC(seq: str) -> float:
    enc = cb.scores.EffectiveNumberOfCodons(bg_correction=True)
    enc_score = enc.get_score(seq) / 61
    return enc_score


def seq_entropy(seq: str) -> float:
    freqs = [seq.count(base) / len(seq) for base in "ACGT"]
    return entropy(freqs) / 2


def count_g_runs(seq: str, min_run_length: int = 4) -> float:
    """
    Calculates the fraction of the sequence that contains G-runs
    of at least 'min_run_length' (like 'GGGG').
    Normalized by the sequence length.
    """
    if len(seq) == 0:
        return 0.0
    seq = seq.upper()
    count = 0
    i = 0
    while i < len(seq):
        if seq[i] == 'G':
            run_length = 1
            while i + 1 < len(seq) and seq[i + 1] == 'G':
                run_length += 1
                i += 1
            if run_length >= min_run_length:
                count += 1
        i += 1
    return count / len(seq)


def hairpin_score(seq: str, min_overlap: int = 4) -> float:
    """
    Estimates the potential of the sequence to form a hairpin structure
    by checking how many small subsequences appear in its reverse complement.
    """
    seq = seq.upper()
    antisense = get_antisense(seq)
    matches = 0
    for i in range(len(seq) - min_overlap + 1):
        sub = seq[i:i + min_overlap]
        if sub in antisense:
            matches += 1
    return matches / len(seq)


def gc_content_3prime_end(aso_sequence: str, window: int = 5) -> float:
    """Calculate the GC content at the 3' end of the ASO sequence."""
    if len(aso_sequence) < window:
        return 0.0
    three_prime_end = aso_sequence[-window:]
    gc_count = three_prime_end.count('G') + three_prime_end.count('C')
    return gc_count / window


def toxic_motif_count(aso_sequence, motifs=['UGU', 'GGTGG', 'TGGT', 'GGGU']) -> float:
    """
    Counts the number of toxic motif appearances in the ASO.
    Returns normalized count (0–1) based on max possible motif hits.

    Parameters:
        aso_sequence (str): DNA or RNA ASO sequence
        motifs (list of str): known toxic motifs (excluding 'GGGG' to avoid overlap with G-runs)

    Returns:
        float: normalized toxic motif count
    """
    sequence = aso_sequence.upper()
    total = 0
    for motif in motifs:
        total += len(re.findall(motif, sequence))

    max_possible = len(sequence)  # conservative upper bound
    return min(total / max_possible, 1.0)


def nucleotide_diversity(seq: str) -> float:
    # checking the nucleotide diversity of the ASO sequence and normalize it by the
    # max value 16
    nucs = [seq[i:i+2] for i in range(len(seq)-1)]
    unique = set(nucs)
    return len(unique) / 16

def stop_codon_count(seq: str, codons=('TAA', 'TAG', 'TGA')) -> float:
    """
    Counts occurrences of stop codons (TAA, TAG, TGA) in the sequence.
    Returns a normalized value (count per length).
    Returns:
        float: normalized count of stop codons
    """
    seq = seq.upper()
    count = sum(seq.count(codon) for codon in codons)
    return count / len(seq)

def tandem_repeats_score(seq: str, min_unit=2, max_unit=6) -> float:
    """
    Calculates how many short motifs (2–6 nt) repeat consecutively (tandem repeats).
    Useful for detecting repetitive structures that may form hairpins or reduce specificity.
    """
    score = 0
    for unit_len in range(min_unit, max_unit + 1):
        for i in range(len(seq) - unit_len * 2 + 1):
            unit = seq[i:i+unit_len]
            repeat_count = 1
            j = i + unit_len
            while j + unit_len <= len(seq) and seq[j:j+unit_len] == unit:
                repeat_count += 1
                j += unit_len
            if repeat_count >= 2:
                score += repeat_count - 1
    return score/len(seq)

def flexible_dinucleotide_fraction(seq: str) -> float:
    """
    Computes the fraction of flexible dinucleotides (AT and TA) in the sequence.
    These dinucleotides are less stable and can indicate structurally flexible regions.
    Args:
        seq (str): DNA sequence (assumed to be uppercase A/C/G/T)
    Returns:
        float: Fraction of AT or TA dinucleotides, normalized by sequence length
    """
    if len(seq) < 2:
        return 0.0
    seq = seq.upper()
    count = 0
    for i in range(len(seq)-1):
        pair = seq[i:i+2]
        if pair in ['AT','TA']:
            count += 1
    return count/(len(seq)-1)

####################################################################
def calculate_chimera_ars(suffix_array, target_sequence, step_size):
    longest_prefix_lengths = []

    for start_index in range(1, len(target_sequence), step_size):
        prefix, _ = longest_prefix(target_sequence[start_index:], suffix_array)
        longest_prefix_lengths.append(len(prefix))

    chimera_ars_score = np.mean(longest_prefix_lengths)
    return chimera_ars_score
###################################################################
