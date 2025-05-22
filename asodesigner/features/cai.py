def calc_CAI_weight(reference_seq):
    # keys =  ['phe','ile','val','pro','thr','ala','Tyr','His','Gln','Asn','Lys','Asp','Glu','Cys','Gly']
    phe_dict = {"TTT": 0, "TTC": 0}
    ile_dict = {"ATT": 0, "ATC": 0, "ATA": 0}
    val_dict = {"GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0}
    pro_dict = {"CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0}
    thr_dict = {"ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0}
    ala_dict = {"GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0}
    tyr_dict = {"TAT": 0, "TAC": 0}
    his_dict = {"CAT": 0, "CAC": 0}
    gln_dict = {"CAA": 0, "CAG": 0}
    asn_dict = {"AAT": 0, "AAC": 0}
    lys_dict = {"AAA": 0, "AAG": 0}
    asp_dict = {"GAT": 0, "GAC": 0}
    glu_dict = {"GAA": 0, "GAG": 0}
    cys_dict = {"TGT": 0, "TGC": 0}
    gly_dict = {"GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}
    leu_dict = {"TTA": 0, "TTG": 0, "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0}
    ser_dict = {"TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0, "AGT": 0, "AGC": 0}
    arg_dict = {"CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0, "AGA": 0, "AGG": 0}
    dict_lst = [phe_dict, ile_dict, val_dict, pro_dict, thr_dict, ala_dict, tyr_dict, his_dict, gln_dict, asn_dict,
                lys_dict, asp_dict, glu_dict, cys_dict, gly_dict, leu_dict, ser_dict, arg_dict]
    # TODO: check if reference_seq is right
    index_lst = list(range(0, len(reference_seq), 3))
    for n in index_lst:
        curr_codon = reference_seq[n:n + 3]
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


def calc_CAI(seq, CAI_weights):
    CAI = 1
    # clean_seq = seq[3:-3] #removing the start and stop codon
    codon_num = len(seq) // 3
    index_lst = list(range(0, len(seq), 3))
    for n in index_lst:
        curr_codon = seq[n:n + 3]
        for dic in CAI_weights:
            for key in dic:
                if key == curr_codon and dic[key] != 0:
                    CAI *= (dic[key] ** (1 / codon_num))
                    # print(CAI)
    return CAI
