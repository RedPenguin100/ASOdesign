from collections import Counter
import pandas as pd
import primer3
from primer3 import calc_hairpin

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

###########################################################################################


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
