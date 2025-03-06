from Bio import SeqIO

def get_gfp_first_exp():
    gfp_obj = next(SeqIO.parse('../data/GFP_first_exp.fasta', 'fasta'))
    gfp_seq = str(gfp_obj.seq.upper())

    with open('../data/GFP_context.txt', 'r') as f:
        gfp_context = f.read().upper()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    gap = 100
    gfp_ext = gfp_context[gfp_start - gap: gfp_start + len(gfp_seq) + gap]

    return gfp_ext