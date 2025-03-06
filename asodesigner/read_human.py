import bisect
import gzip
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
from asodesigner.consts import HUMAN_TRANSCRIPTS_FASTA_GZ, HUMAN_GTF_BASIC, HUMAN_GFF, HUMAN_DB_BASIC, \
    HUMAN_GENOME_FASTA_GZ, \
    HUMAN_GENOME_FASTA, HUMAN_DB_BASIC_INTRONS, HUMAN_TRANSCRIPTS_FASTA


def cond_print(text, verbose=False):
    if verbose:
        print(text)


class LocusInfo:
    def __init__(self):
        self.exons = []
        self.introns = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""


def get_locus_to_data_dict():
    db_path = Path(HUMAN_DB_BASIC_INTRONS)

    if not db_path.exists():
        start = time.time()
        db = gffutils.create_db(str(HUMAN_GFF), dbfn=str(db_path), force=True, keep_order=True,
                                merge_strategy='merge', sort_attribute_values=True)
        end = time.time()
        print(f"DB create took: {end - start}s")
    else:
        db = gffutils.FeatureDB(str(db_path))
        # db.update(list(db.create_introns())) # Uncomment to create introns

    start = time.time()
    # with gzip.open(str(HUMAN_GENOME_FASTA_GZ), 'rt') as handle:
    #     fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    with open(str(HUMAN_TRANSCRIPTS_FASTA), 'r') as handle:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    end = time.time()
    print(f"Time to unzip fasta: {end - start}s")

    locus_to_data = dict()
    locus_to_strand = dict()

    for feature in db.features_of_type(('CDS', 'UTR'), order_by='start'):
        chrom = feature.seqid
        if 'chrM' == chrom:
            continue
        locus_tags = feature.attributes['gene_name']
        if len(locus_tags) != 1:
            raise ValueError(f"Multiple locuses: {locus_tags}")
        locus_tag = locus_tags[0]

        if feature.featuretype == 'CDS':
            cds = feature
            seq = fasta_dict[chrom].seq[cds.start - 1: cds.end]
            if cds.strand == '-':
                seq = seq.reverse_complement()
            seq = seq.upper()

            if locus_tag not in locus_to_data:
                locus_info = LocusInfo()
                locus_info.exons = [(cds.start, seq)]
                locus_to_data[locus_tag] = locus_info

                locus_to_strand[locus_tag] = cds.strand
            else:
                bisect.insort(locus_to_data[locus_tag].exons, (cds.start, seq))

    for locus_tag in locus_to_data:
        l = locus_to_data[locus_tag].exons
        if locus_to_strand[locus_tag] == '-':
            l.reverse()

        locus_to_data[locus_tag].exons = [element for _, element in l]
    print(f"Reversing took: {end - start}s")

    return locus_to_data


if __name__ == '__main__':
    start = time.time()
    gene_to_data = get_locus_to_data_dict()
    end = time.time()
    print(f"Time to read full human: {end - start}s")

    start = time.time()
    i = 0
    for gene in gene_to_data.items():
        print(gene[0])
        i += len(gene[1].exons[0])
        continue
    end = time.time()
    print(f"Iterate took: {end - start}s, i={i}")
    # print(gene_to_data)
    # print(len(gene_to_data))
