import time

from asodesigner.read_yeast import get_locus_to_data_dict_alternative, get_locus_to_data_dict
import difflib
import pytest


def test_locus():
    start = time.time()
    locus_to_data = get_locus_to_data_dict()
    end = time.time()
    print(f"Regular Took: {end - start}s")

    start = time.time()
    locus_to_data_alt = get_locus_to_data_dict_alternative()
    end = time.time()
    print(f"DB Took: {end - start}s")

    for key, value in locus_to_data.items():
        alt_exons = locus_to_data_alt[key].exons
        exons = locus_to_data[key].exons
        assert len(exons) == len(alt_exons)

        for i in range(len(exons)):
            exon = exons[i]
            alt_exon = alt_exons[i]

            assert exon == alt_exon, "key: " + key + " len exons " + str(len(exons)) + "i: " + str(i)


def test_intron_regression():
    locus_to_data = get_locus_to_data_dict_alternative()

    assert len(locus_to_data['YNCA0002W'].introns) == 1
    intron = locus_to_data['YNCA0002W'].introns[0]
    assert str(intron) == 'CGACTTCCTGATTAAACAGGAAGACAAAGCA'