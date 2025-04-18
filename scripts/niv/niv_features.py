from collections import Counter
import pandas as pd
import codonbias as cb
import primer3
from primer3 import calc_hairpin
from suffix_array import calc_suffix_array, longest_prefix
import numpy as np
from scipy.stats import entropy
import math as math
import re
from Bio.Seq import Seq
################################################################
def compute_ENC(seq):
    enc = cb.scores.EffectiveNumberOfCodons(bg_correction=True)
    enc_score = enc.get_score(seq) / 61
    return enc_score

##################################################################
def palindromic_fraction(seq, l):
    def is_palindrome(seq):
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

    count = 0
    if ~l % 2:
        for n in range(len(seq) - l + 1):
            curr_seq = seq[n:n + l]
            count += is_palindrome(curr_seq)
    return count / len(seq)

def homooligo_count (seq):
    seq +='$'
    tot_count = 0
    curr_seq = ''
    n = 0
    while n in range(len(seq)-1):
        while seq[n] == seq[n+1]:
            curr_seq += seq[n]
            n += 1
        if len(curr_seq) > 1:
            tot_count += len(curr_seq)+1
        n += 1
        curr_seq = ''
    return tot_count/len(seq)
####################################################################
def calculate_chimera_ars(suffix_array, target_sequence, step_size):
        longest_prefix_lengths = []

        for start_index in range(1, len(target_sequence), step_size):
            prefix, _ = longest_prefix(target_sequence[start_index:], suffix_array)
            longest_prefix_lengths.append(len(prefix))

        chimera_ars_score = np.mean(longest_prefix_lengths)
        return chimera_ars_score
    ####################
###################################################################
def calc_CAI_weight(reference_seq):
    #keys =  ['phe','ile','val','pro','thr','ala','Tyr','His','Gln','Asn','Lys','Asp','Glu','Cys','Gly']
    phe_dict =  {"TTT": 0, "TTC":0}
    ile_dict = {"ATT":0, "ATC":0, "ATA":0}
    val_dict = {"GTT":0, "GTC":0, "GTA":0, "GTG":0}
    pro_dict = {"CCT":0, "CCC":0, "CCA":0, "CCG":0}
    thr_dict = {"ACT":0, "ACC":0, "ACA":0, "ACG":0}
    ala_dict = {"GCT":0, "GCC":0, "GCA":0, "GCG":0}
    tyr_dict = {"TAT":0, "TAC":0}
    his_dict = {"CAT":0, "CAC":0}
    gln_dict = {"CAA":0, "CAG":0}
    asn_dict = {"AAT":0, "AAC":0}
    lys_dict = {"AAA":0, "AAG":0}
    asp_dict = {"GAT":0, "GAC":0}
    glu_dict = {"GAA":0, "GAG":0}
    cys_dict = {"TGT":0, "TGC":0}
    gly_dict = {"GGT":0, "GGC":0, "GGA":0, "GGG":0}
    leu_dict = {"TTA":0, "TTG":0, "CTT":0, "CTC":0, "CTA":0, "CTG":0}
    ser_dict = {"TCT":0, "TCC":0, "TCA":0, "TCG":0, "AGT":0, "AGC":0}
    arg_dict = {"CGT":0, "CGC":0, "CGA":0, "CGG":0, "AGA":0, "AGG":0}
    dict_lst = [phe_dict, ile_dict, val_dict, pro_dict, thr_dict, ala_dict, tyr_dict, his_dict, gln_dict, asn_dict, lys_dict, asp_dict, glu_dict, cys_dict, gly_dict,leu_dict,ser_dict,arg_dict]
    index_lst = list(range(0,len(seq),3))
    for n in index_lst:
        curr_codon = reference_seq[n:n+3]
        for dic in dict_lst:
            for key in dic:
                if curr_codon == key:
                    dic[key] += 1
    w_dict_lst = []
    for dic in dict_lst:
        max_val = max(dic.values())
        if max_val != 0:
            dic = {key: value / max_val for key, value in dic.items()}
        w_dict_lst.append(dic)
    return w_dict_lst
########################################################################
def calc_CAI(seq, CAI_weights):
    CAI = 1
    #clean_seq = seq[3:-3] #removing the start and stop codon
    codon_num = len(seq)//3
    index_lst = list(range(0,len(seq),3))
    for n in index_lst:
        curr_codon = seq[n:n+3]
        for dic in CAI_weights:
            for key in dic:
                if key == curr_codon and dic[key]!=0:
                    CAI *= (dic[key]** (1/codon_num))
                    #print(CAI)
    return CAI

################################################################
def seq_entropy(seq):
    freqs = [seq.count(base) / len(seq) for base in "ACGT"]
    return entropy(freqs)/2

###############################################################
def calc_tAI(seq, weight_dictionary):
    index_lst = list(range(0,len(seq),3))
    tAI_log = 0
    L = len(index_lst)
    for n in index_lst:
        curr_codon = seq[n:n+3]
        if len(curr_codon)!= 3:
            continue
        else:
            tAI = weight_dictionary[curr_codon]
        if tAI != 0:
           tAI_log += math.log(tAI)
    return math.exp(tAI_log/L)
#############################################################
def tai_weights(category):
    # each category will give you the weights according to the creature you chose
    # please follow the next rules:
    # sc: will give yot the weight for Saccharomyces cerevisiae

    if category == "sc":
        weight_dict = {}
        weight_dict["AAA"] = 7; weight_dict["AAC"] = 10; weight_dict["AAG"] = 16.24; weight_dict["AAT"] = 5.9
        weight_dict["ACA"] = 4.0011; weight_dict["ACC"] = 7.92; weight_dict["ACG"] = 2.28; weight_dict["ACT"] = 11
        weight_dict["AGA"] = 11; weight_dict["AGC"] = 4; weight_dict["AGG"] = 4.52; weight_dict["AGT"] = 2.36
        weight_dict["ATA"] = 2.0013; weight_dict["ATC"] = 9.36; weight_dict["ATG"] = 10.64; weight_dict["ATT"] = 18.9

        weight_dict["CAA"] = 9; weight_dict["CAC"] = 7; weight_dict["CAG"] = 3.88; weight_dict["CAT"] = 4.13
        weight_dict["CCA"] = 10.0002; weight_dict["CCC"] = 1.44; weight_dict["CCG"] = 3.2; weight_dict["CCT"] = 2
        weight_dict["CGA"] = 0.0006; weight_dict["CGC"] = 4.32; weight_dict["CGG"] = 1; weight_dict["CGT"] = 6
        weight_dict["CTA"] = 3; weight_dict["CTC"] = 1; weight_dict["CTG"] = 0.96; weight_dict["CTT"] = 0.59

        weight_dict["GAA"] = 14; weight_dict["GAC"] = 16; weight_dict["GAG"] = 6.48; weight_dict["GAT"] = 9.44
        weight_dict["GCA"] = 5.0011; weight_dict["GCC"] = 7.92; weight_dict["GCG"] = 1.6; weight_dict["GCT"] = 11
        weight_dict["GGA"] = 3; weight_dict["GGC"] = 16; weight_dict["GGG"] = 2.96; weight_dict["GGT"] = 9.44
        weight_dict["GTA"] = 2.0014; weight_dict["GTC"] = 10.08; weight_dict["GTG"] = 2.64; weight_dict["GTT"] = 14

        weight_dict["TAA"] = 0; weight_dict["TAC"] = 8; weight_dict["TAG"] = 0; weight_dict["TAT"] = 4.72
        weight_dict["TCA"] = 3.0011; weight_dict["TCC"] = 7.92; weight_dict["TCG"] = 1.96; weight_dict["TCT"] = 11
        weight_dict["TGA"] = 0; weight_dict["TGC"] = 4; weight_dict["TGG"] = 6; weight_dict["TGT"] = 2.36
        weight_dict["TTA"] = 7; weight_dict["TTC"] = 10; weight_dict["TTG"] = 12.24; weight_dict["TTT"] = 5.9
        weight_dict = {k: v / max(weight_dict.values()) for k, v in weight_dict.items()}
        return  weight_dict
####################################################################################################################################

def purine_content(seq):
    """
    Calculates the fraction of purine bases (A and G) in the sequence.
    Purine-rich sequences may be more stable and bind better to RNA targets.
    """
    seq = seq.upper()
    count = seq.count("A") + seq.count("G")
    if len(seq) == 0:
        return 0.0
    return count / len(seq)
#################################################################

def count_g_runs(seq, min_run_length=4):
    """
    Calculates the fraction of the sequence that contains G-runs
    of at least 'min_run_length' (like 'GGGG').
    Normalized by the sequence length.
    """
    if len(seq) == 0:
        return  0.0
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
#############################################################

def hairpin_score(seq: str, min_overlap: int = 4):
    """
    Estimates the potential of the sequence to form a hairpin structure
    by checking how many small subsequences appear in its reverse complement.
    """
    seq = seq.upper()
    rc = str(Seq(seq).reverse_complement())
    matches = 0
    for i in range(len(seq) - min_overlap + 1):
        sub = seq[i:i+min_overlap]
        if sub in rc:
            matches += 1
    return matches/len(seq)
########################################################################################
def gc_content_3prime_end(aso_sequence, window=5):
    """Calculate the GC content at the 3' end of the ASO sequence."""
    if len(aso_sequence) < window:
        return 0.0
    three_prime_end = aso_sequence[-window:]
    gc_count = three_prime_end.count('G') + three_prime_end.count('C')
    return gc_count/window
########################################################################################

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
#######################################################################3

def hairpin_dG_energy(seq: str):
    """
    Returns the raw ΔG (Gibbs free energy) of predicted hairpin structure,
    divided by sequence length.

    This value is NOT normalized to [0,1].
    Positive values suggest unstable/no structure.
    Negative values indicate stronger/stable hairpins that may interfere with ASO activity.
    """
    hairpin = calc_hairpin(seq)
    print("structure_found:", hairpin.structure_found)  # ← נוספה שורה לבדיקת הדגל
    if not hairpin.structure_found:
        return 0
    return hairpin.dg if len(seq) > 0 else 0


################################################################

def hairpin_tm(seq: str):
    """
    Feature: hairpin_tm
    Melting temperature (Tm) of the predicted hairpin structure.
    Higher Tm means the structure is more stable at physiological temperatures.

    Returns:
        float: Tm of the hairpin
    """
    if len(seq) == 0:
        return 0
    hairpin = primer3.calc_hairpin(seq)
    if not hairpin.structure_found:
        return 0
    return hairpin.tm

#########################################################
def tandem_repeats_score(seq, min_unit=2, max_unit=6):
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

################################################################
def dispersed_repeats_score(seq, min_unit=2, max_unit=6):
    """
    Counts motifs (2–6 nt) that appear more than once, even if not consecutive.
    Helps detect internal similarity and potential self-binding regions.
    """
    unit_counter = Counter()
    for unit_len in range(min_unit, max_unit + 1):
        for i in range(len(seq) - unit_len + 1):
            unit = seq[i:i+unit_len]
            unit_counter[unit] += 1
    score = sum(count - 1 for count in unit_counter.values() if count > 1)
    return score / len(seq)
#################################################################
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
########################################################################
def gc_skew(seq: str) -> float:
    """
    Computes GC skew = (G - C) / (G + C)
    A measure of strand asymmetry in G/C content, can affect hybridization and folding.
    """
    seq = seq.upper()
    G_counts = seq.count("G")
    C_counts = seq.count("C")
    if G_counts + C_counts == 0:
        return  0.0
    return (G_counts - C_counts) / (G_counts + C_counts)
#############################################################################################
def gc_skew_ends(seq: str, window: int = 5) -> float:
    """
    Calculates the GC-content difference between the 5' and 3' ends of the sequence.
    Measures thermodynamic asymmetry between ends.
    """
    seq = seq.upper()
    if len(seq) < 2 * window:
        return 0.0
    start = seq[:window]
    end = seq[-window:]
    gc_5 = start.count("G") + start.count("C")
    gc_3 = end.count("G") + end.count("C")
    return (gc_5 - gc_3) / window
##################################################################
def at_skew(seq: str) -> float:
    """
    Calculates AT skew = (A - T) / (A + T)
    A measure of asymmetry in A/T content, can affect flexibility and binding dynamics.
    """
    seq = seq.upper()
    A_counts = seq.count("A")
    T_counts = seq.count("T")
    if A_counts + T_counts == 0:
        return 0.0  # avoid division by zero
    return (A_counts - T_counts) / (A_counts + T_counts)
####################################################################
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
###########################################################################################

def nucleotide_diversity(seq: str) -> float:
    # checking the nucleotide diversity of the ASO sequence and normalize it by the
    # max value 16
    nucs = [seq[i:i+2] for i in range(len(seq)-1)]
    unique = set(nucs)
    return len(unique) / 16
###################################################################################

df = pd.read_csv("ANGPTL2_antisense.csv")
dict_tai = tai_weights("sc")
df['tai']= df["Sequence"].apply(lambda x:calc_tAI(x,dict_tai))
df["purine content"]=df["Sequence"].apply(lambda x: purine_content(x))
df["gggg content"] = df['Sequence'].apply(lambda x: count_g_runs(x))
df["hairpin score"] = df['Sequence'].apply(lambda x: hairpin_score(x))
df["gc 3 prime"] = df["Sequence"].apply(lambda x:gc_content_3prime_end(x) )
df["toxic motif"] = df["Sequence"].apply(lambda x:toxic_motif_count(x))
df[" tandem repeats"] = df["Sequence"].apply(lambda  x: tandem_repeats_score(x))
df["flexible_dinucleotide_fraction"] = df["Sequence"].apply(lambda x: flexible_dinucleotide_fraction(x))
df["repeats"] = df["Sequence"]. apply(lambda x : dispersed_repeats_score(x))
df["gc_skew"] = df["Sequence"].apply(lambda x : gc_skew(x))
df["gc_skew_ends"] =df["Sequence"].apply(lambda x:gc_skew_ends(x))
print(df)

seq = "ATGCGTACGTAGCTAGCTA"

hairpin = primer3.calc_hairpin(seq)
print("Hairpin ΔG normalized:", hairpin_dG_energy(seq))
print("Hairpin Tm:", hairpin_tm(seq))
