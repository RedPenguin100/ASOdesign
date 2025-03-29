import bisect
from pathlib import Path

import gffutils
from numba import njit
from numba.typed import Dict

from asodesigner.consts import HUMAN_GFF, HUMAN_DB_BASIC_INTRONS, HUMAN_DB_BASIC_INTRONS_GZ
from asodesigner.file_utils import read_human_genome_fasta_dict
from asodesigner.process_utils import LocusInfo
from asodesigner.timer import Timer


def cond_print(text, verbose=False):
    if verbose:
        print(text)

def create_human_genome_db(path: Path, create_introns=False):
    print("Creating human genome database. WARNING - this is slow!")
    with Timer() as t:
        db = gffutils.create_db(str(HUMAN_GFF), dbfn=str(path), force=True, keep_order=True,
                                merge_strategy='merge', sort_attribute_values=True)
        if create_introns:
            db.update(list(db.create_introns()))
    print(f"DB create took: {t.elapsed_time}s")
    return db


def get_human_genome_annotation_db(create_db=False):
    db_path = HUMAN_DB_BASIC_INTRONS

    if not db_path.is_file():
        if HUMAN_DB_BASIC_INTRONS_GZ.is_file():
            raise ValueError(
                f"DB file is not unzipped: {HUMAN_DB_BASIC_INTRONS_GZ}, please unzip to use! (Consider README.md)")

        if create_db:
            db = create_human_genome_db(db_path, create_introns=True)
        else:
            raise ValueError(
                f"DB not found in path: {str(db_path)}, either download it or create (please consider README.md)")
    else:
        db = gffutils.FeatureDB(str(db_path))
    return db


def get_locus_to_data_dict(create_db=False, include_introns=False):
    db = get_human_genome_annotation_db(create_db)
    fasta_dict = read_human_genome_fasta_dict()

    locus_to_data = dict()
    locus_to_strand = dict()

    if include_introns:
        feature_types = ('exon', 'intron')
    else:
        feature_types = ('exon')


    for feature in db.features_of_type(feature_types, order_by='start'):
        chrom = feature.seqid
        if 'chrM' == chrom:
            continue
        locus_tags = feature.attributes['gene_name']
        if len(locus_tags) != 1:
            raise ValueError(f"Multiple loci: {locus_tags}")
        locus_tag = locus_tags[0]

        if feature.featuretype == 'exon':
            cds = feature
            seq = fasta_dict[chrom].seq[cds.start - 1: cds.end]
            if cds.strand == '-':
                seq = seq.reverse_complement()
            seq = seq.upper()

            if locus_tag not in locus_to_data:
                locus_info = LocusInfo()
                locus_info.exons = [(cds.start, seq)]
                locus_info.introns = []
                locus_to_data[locus_tag] = locus_info

                locus_to_strand[locus_tag] = cds.strand
            else:
                bisect.insort(locus_to_data[locus_tag].exons, (cds.start, seq))
        elif feature.featuretype == 'intron' and include_introns:
            seq = fasta_dict[chrom].seq[feature.start - 1: feature.end]

            if feature.strand == '-':
                seq = seq.reverse_complement()
            seq = seq.upper()

            intron = feature

            if locus_tag not in locus_to_data:
                locus_info = LocusInfo()
                locus_info.introns = [(feature.start, seq)]
                locus_info.exons = []
                locus_to_data[locus_tag] = locus_info

                locus_to_strand[locus_tag] = intron.strand
            else:
                bisect.insort(locus_to_data[locus_tag].introns, (intron.start, seq))
        elif feature.featuretype == 'UTR':
            pass
        else:
            print("Feature type: ", feature.featuretype)


    for locus_tag in locus_to_data:
        locus_info = locus_to_data[locus_tag]
        if locus_to_strand[locus_tag] == '-':
            locus_info.exons.reverse()
            if include_introns:
                locus_info.introns.reverse()
        locus_info.exons = [element for _, element in locus_info.exons]

        if include_introns:
            locus_info.introns = [element for _, element in locus_info.introns]

    return locus_to_data


if __name__ == '__main__':
    with Timer() as t:
        gene_to_data = get_locus_to_data_dict(include_introns=True)
    print(f"Time to read full human: {t.elapsed_time}s")

    with Timer() as t:
        i = 0
        for gene in gene_to_data.items():
            if len(gene[1].exons) == 0:
                print("Weird gene, ", gene[0])
            else:
                i += len(gene[1].exons[0])
            continue
    print(f"Iterate took: {t.elapsed_time}s, i={i}")
    # print(gene_to_data)
    # print(len(gene_to_data))
