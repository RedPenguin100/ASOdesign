import math


def calc_tAI(seq : str, weight_dictionary) -> float:
    index_lst = list(range(0, len(seq), 3))
    tAI_log = 0
    L = len(index_lst)
    for n in index_lst:
        curr_codon = seq[n:n + 3]
        curr_codon = curr_codon.replace('U', 'T')
        if len(curr_codon) != 3:
            continue
        else:
            tAI = weight_dictionary[curr_codon]
        if tAI != 0:
            tAI_log += math.log(tAI)
    return math.exp(tAI_log / L)


def tai_weights(category):
    # each category will give you the weights according to the creature you chose
    # please follow the next rules:
    # sc: will give yot the weight for Saccharomyces cerevisiae

    if category == "sc":
        weight_dict = dict()

        # @formatter:off
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
        # @formatter:on

        weight_dict = {k: v / max(weight_dict.values()) for k, v in weight_dict.items()}

        return weight_dict
