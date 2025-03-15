import bisect
import gzip
from pathlib import Path

import gffutils
from Bio import SeqIO
from asodesigner.consts import HUMAN_GFF, HUMAN_DB_BASIC_INTRONS, HUMAN_GENOME_FASTA
from asodesigner.file_utils import get_fasta_dict_from_path
from asodesigner.timer import Timer


def cond_print(text, verbose=False):
    if verbose:
        print(text)


class LocusInfo:
    def __init__(self):
        self.exons = []
        self.introns = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""


def create_human_genome_db(path: Path):
    print("Creating human genome database. WARNING - this is slow!")
    with Timer() as t:
        db = gffutils.create_db(str(HUMAN_GFF), dbfn=str(path), force=True, keep_order=True,
                                merge_strategy='merge', sort_attribute_values=True)
    print(f"DB create took: {t.elapsed_time}s")
    return db


def get_locus_to_data_dict():
    db_path = Path(HUMAN_DB_BASIC_INTRONS)

    if not db_path.exists():
        # db = create_human_genome_db(db_path)
        raise ValueError(f"Please read the README.md! After that put the db path in {str(db_path)}")
    else:
        db = gffutils.FeatureDB(str(db_path))
        # db.update(list(db.create_introns())) # Uncomment to create introns

    with Timer() as t:
        fasta_dict = get_fasta_dict_from_path(HUMAN_GENOME_FASTA)

    print(f"Time to unzip fasta: {t.elapsed_time}s")

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

    return locus_to_data


if __name__ == '__main__':
    with Timer() as t:
        gene_to_data = get_locus_to_data_dict()

    print(f"Time to read full human: {t.elapsed_time}s")

    with Timer() as t:
        i = 0
        for gene in gene_to_data.items():
            print(gene[0])
            i += len(gene[1].exons[0])
            continue
    print(f"Iterate took: {t.elapsed_time}s, i={i}")
    # print(gene_to_data)
    # print(len(gene_to_data))
