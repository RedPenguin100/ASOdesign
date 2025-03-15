import pytest

from asodesigner.fold import get_mfe_scores, get_trigger_mfe_scores_by_risearch, Interaction, dump_target_file
from asodesigner.target_finder import get_gfp_first_exp
from asodesigner.util import get_antisense


def test_risearch():
    # example
    name_to_sequence = {
        "0": "AUGUGUUCAAUUUUUGGCGUAUUCGAUAUCAAAACAGACGCAGUUGAGCUGCGUAAGAAAGCACUCGAGCUGUCACGCCUGAUGCGUCAUCGUGGCCCGGACUGGUCCGGUAUUUAUGCCAGCGAUAACGCCAUUCUCGCCCACGAACGUCUGUCAAUUGUUGACGUUAACGCGGGGGCGCAACCUCUCUACAACCAACAAAAAACUCAUGUGCUGGCGGUAAACGGUGAAAUCUACAACCACCAGGCACUGCGCGCCGAAUAUGGCGAUCGUUAUCAGUUCCAGACCGGAUCUGACUGUGAAGUGAUCCUCGCGCUGUAUCAGGAAAAAGGGCCGGAAUUUCUCGACGACUUGCAGGGCAUGUUUGCCUUUGCCUUGUACGACAGCGAAAAAGAUGCCUACCUGAUUGGUCGCGACCAUCUGGGGAUCAUCCCACUGUAUAUGGGCUAUGACGAACACGGUCAGCUGUAUGUGGCCUCAGAAAUGAAAGCCCUGGUGCCAGUUUGCCGCACGAUUAAAGAGUUCCCGGCGGGGAGCUAUUUGUGGAGCCAGGACGGCGAAAUCCGUUCUUAUUAUCAUCGCGACUGGUUCGACUACGAUGCGGUGAAAGAUAACGUAACCGACAAAAACGAGCUGCGUCAGGCACUUGAAGAUUCCGUUAAAAGCCAUCUGAUGUCUGAUGUGCCUUACGGUGUGCUGCUUUCUGGUGGUCUGGAUUCCUCAAUUAUUUCCGCUAUCACCAAGAAAUACGCAGCCCGUCGCGUGGAAGAUCAGGAACGCUCUGAAGCCUGGUGGCCGCAGUUACACUCCUUUGCUGUAGGUCUGCCGGGUUCACCGGAUCUUAAAGCAGCCCAGGAAGUGGCAAACCAUCUGGGCACGGUGCAUCACGAAAUUCACUUCACUGUACAGGAAGGUCUGGAUGCCAUCCGCGACGUGAUUUACCACAUCGAAACUUAUGAUGUGACCACAAUUCGCGCUUCAACACCGAUGUAUUUAAUGUCGCGUAAGAUCAAGGCGAUGGGCAUUAAAAUGGUGCUGUCCGGUGAAGGUUCUGAUGAAGUGUUUGGCGGUUAUCUUUACUUCCAUAAAGCGCCCAACGCUAAAGAACUGCAUGAAGAGACGGUGCGUAAACUGCUGGCCCUGCAUAUGUAUGACUGCGCGCGCGCCAACAAAGCGAUGUCAGCCUGGGGCGUGGAAGCACGCGUUCCGUUCCUCGACAAAAAAUUCCUUGAUGUGGCGAUGCGCAUUAACCCGCAGGAUAAAAUGUGCGGUAACGGCAAAAUGGAAAAACACAUCCUGCGUGAAUGUUUUGAGUCAUACCUGCCCGCAAGCGUGGCCUGGCGGCAGAAAGAGCAGUUCUCCGAUGGCGUCGGUUACAGUUGGAUCGACACCCUGAAAGAAGUGGCGGCGCAGCAGGUUUCUGAUCAGCAACUGGAAACUGCCCGCUUCCGCUUCCCGUACAACACGCCAACCUCAAAAGAAGCGUAUCUGUACCGGGAGAUCUUUGAAGAACUGUUCCCGCUUCCGAGCGCCGCUGAGUGCGUGCCUGGCGGUCCUUCCGUCGCGUGUUCUUCCGCUAAAGCGAUCGAAUGGGAUGAAGCGUUCAAGAAAAUGGACGAUCCAUCUGGUCGUGCGGUUGGUGUUCACCAGUCGGCAUAUAAGUAA",
        "1": "AUGUUCGAACAACGCGUAAAUUCUGACGUACUGACCGUUUCUACCGUUAACUCUCAGGAUCAGGUAACCCAAAAACCCCUGCGUGACUCGGUUAAACAGGCACUGAAGAACUAUUUUGCUCAACUGAAUGGUCAGGAUGUGAAUGACCUCUAUGAGCUGGUACUGGCUGAAGUAGAACAGCCCCUGUUGGACAUGGUGAUGCAAUACACCCGUGGUAACCAGACCCGUGCUGCGCUGAUGAUGGGCAUCAACCGUGGUACGCUGCGUAAAAAAUUGAAAAAAUACGGCAUGAACUAA"
    }

    target_path = dump_target_file('target-cache.fa', name_to_sequence)
    result = get_trigger_mfe_scores_by_risearch("UAGAUGCGCCACUUGUGGUAUUCCCGCAUC", name_to_sequence, minimum_score=900,
                                                parsing_type='2', target_file_cache=str(target_path))
    print(result)
    mfe_scores = get_mfe_scores(result, '2')

    print(mfe_scores)


def test_bad_fit():
    sense = 'TTTTTTTCTTCCATT'

    result = get_trigger_mfe_scores_by_risearch(sense, {'0': sense + sense}, parsing_type='2', minimum_score=900)
    print(result)

    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[0:15]
    result = get_trigger_mfe_scores_by_risearch('T' + sample_seq[1:15], {"gfp": gfp_seq}, parsing_type='2',
                                                minimum_score=900)
    print(result)


def test_risearch_gfp():
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    print("GFP ontarget(?) : ", gfp_seq[695:714])
    print("Sample", sample_seq)
    print("Sample antisense", get_antisense(sample_seq))
    # name_to_seq = {f"gfp_seq{i}" : gfp_seq for i in range(100)}
    name_to_seq = {f"gfp_seq": gfp_seq}
    result = get_trigger_mfe_scores_by_risearch(sample_seq, name_to_seq,
                                                interaction_type=Interaction.DNA_RNA_NO_WOBBLE,
                                                minimum_score=900, neighborhood=30, parsing_type='2')
    print(result)

    mfe_scores = get_mfe_scores(result, '2')
    print(mfe_scores)

    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    for bad_sample in bad_samples:
        result = get_trigger_mfe_scores_by_risearch(bad_sample, name_to_seq,
                                                    interaction_type=Interaction.DNA_RNA_NO_WOBBLE, minimum_score=900,
                                                    neighborhood=30, parsing_type='2')
        print(result)

        mfe_scores = get_mfe_scores(result, '2')
        print(mfe_scores)
