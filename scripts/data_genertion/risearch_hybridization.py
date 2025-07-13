import random

from scripts.data_genertion.consts import *
from asodesigner.util import get_antisense
from asodesigner.fold import dump_target_file, Interaction, get_trigger_mfe_scores_by_risearch, get_mfe_scores
from modified_dsm import make_dsm_dna_rna, experimental_ps_dna_rna, experimental_neg_ps_dna_rna
from itertools import product


class RNA_TABLES:
    SU95_WOGU = 'dsm_su95_rev_woGU_pos'
    SU95_WGU = 'dsm_su95_rev_wGU_pos'
    SLH04_WOGU = 'dsm_slh04_woGU_pos'


RNA_TABLES = [RNA_TABLES.SU95_WGU, RNA_TABLES.SLH04_WOGU, RNA_TABLES.SU95_WOGU]


# all_data_human_gene_premrna = all_data_human_gene_premrna.sample(frac=0.1)

# genes_u.remove('HIF1A') # the same sequence for all ASOs, only differing chemical modifications
# for base in [118, 122, 126, 130, 134]:
# for dsm in ['dsm_su95_rev_woGU_pos', 'dsm_su95_rev_wGU_pos', 'dsm_slh04_woGU_pos']:

def populate_with_risearch_hybridization(df, genes_u, gene_to_data):
    minimal_scores = [600]

    settings = list(product(RNA_TABLES, [138], minimal_scores, [37], [True, False]))
    for dsm, base, score, temp, transpose in settings:
        print(f"Running {dsm} for base {base} and score {score} temp {temp} transpose {transpose}")

        for gene in genes_u:
            target = gene_to_data[gene].full_mrna
            name_to_sequence = {'target_seq': target}
            parsing_type = '2'

            # prepare and dump target cache
            hash = random.getrandbits(64)
            target_cache_filename = f'target-cache-{hash}.fa'
            target_cache_path = dump_target_file(target_cache_filename, name_to_sequence)

            gene_rows = df[df[CANONICAL_GENE] == gene]

            loop_args = [
                (Interaction.MODIFIED, make_dsm_dna_rna, '2'),
                (Interaction.MODIFIED, experimental_ps_dna_rna, '4'),
                (Interaction.MODIFIED, experimental_neg_ps_dna_rna, '6'),
            ]

            for interaction, func, prefix in loop_args:
                # call your setup function
                if prefix == '2':
                    func(dsm_base=dsm)
                elif prefix == '4':
                    func(base=base, dsm_base=dsm, temp=temp, transpose=transpose)
                elif prefix == '6':
                    func(base=base, dsm_base=dsm, temp=temp, transpose=transpose)
                else:
                    func()

                def _score_one(args):
                    index, row = args
                    antisense = row[SEQUENCE]
                    sense = get_antisense(antisense)
                    l = row[SENSE_LENGTH]
                    tmp = get_trigger_mfe_scores_by_risearch(
                        sense, {'target_seq': target},
                        interaction_type=interaction,
                        minimum_score=score,
                        neighborhood=l,
                        parsing_type=parsing_type,
                        target_file_cache=target_cache_path
                    )
                    scores = get_mfe_scores(tmp, parsing_type)
                    if not scores or not scores[0]:
                        return index, 0, 0.0, 0.0
                    ts = scores[0]
                    return index, len(ts), sum(ts), min(ts)

                # **sequential** scoring loop
                for idx, row in gene_rows.iterrows():
                    idx, cnt, sm, mn = _score_one((idx, row))
                    prefix_col = f"{dsm}{base}{prefix}t{temp}{transpose}on_target_energy_"
                    df.at[idx, prefix_col + f'fits{score}'] = cnt
                    df.at[idx, prefix_col + f'sum{score}']  = sm
                    df.at[idx, prefix_col + f'max{score}']  = mn
