import pytest
import pandas as pd
from Bio.Seq import Seq

from alternative_splicing import (
    get_transcripts_of_gene, get_sub_df, find_in_transcripts,
    find_aso, get_id_from_annotations, get_seq_for_genes
)

# ---------- MOCKING HELPERS ----------
class DummyGene:
    def __init__(self, seq):
        self.full_mrna = seq

# ---------- TESTS ----------
def test_get_transcripts_of_gene(monkeypatch):
    def mock_seq(gene_id, isoforms):
        return ['>ENST000001', 'AUGCUACGGAU', '>ENST000002', 'GGCUACAGUUU']
    monkeypatch.setattr("gget.seq", mock_seq)
    result = get_transcripts_of_gene("ENSG000001")
    assert isinstance(result, list)
    assert len(result) == 4
    assert result[0] == '>ENST000001'

def test_get_sub_df():
    df = pd.DataFrame({
        "Sequence": ["AAA", "AAA", "CCC", "GGG"],
        "Canonical Gene Name": ["GENE1", "GENE1", "GENE2", "GENE3"]
    })
    gene_list = ["GENE1", "GENE3"]
    result = get_sub_df(df, gene_list)
    assert len(result) == 2
    assert "Sequence" in result.columns
    assert "Canonical Gene Name" in result.columns

def test_find_in_transcripts(monkeypatch):
    def mock_find_if_seq_in_gene(seq, transcript):
        return transcript.find(seq)

    monkeypatch.setattr("scripts.data_genertion.extracting_data_from_article.find_if_seq_in_gene", mock_find_if_seq_in_gene)

    transcripts_dict = {
        "GENE1": [">ENST000001", "CCCAAAGGGTTT", ">ENST000002", "TTTAAAGGGCCC"]
    }
    pos, transcript_id, transcript_seq = find_in_transcripts("AAAGGG", "GENE1", transcripts_dict)
    assert pos == 3
    assert transcript_id == "ENST000001"

def test_find_aso(monkeypatch):
    df = pd.DataFrame({
        "Sequence": ["TTT"],
        "Canonical Gene Name": ["GENE1"]
    })
    locus_to_data = {"GENE1": "dummy_id"}
    dict_with_seq = {"GENE1": "AAATTTGGG"}

    def mock_find(seq, full_seq):
        return full_seq.find(seq)

    def mock_transcripts(_):
        return [">TX1", "GGGTTTCCC"]

    monkeypatch.setattr("scripts.data_genertion.extracting_data_from_article.find_if_seq_in_gene", mock_find)
    monkeypatch.setattr("gget.seq", lambda x, isoforms: mock_transcripts(x))

    df_out, transcript_dict = find_aso(df, locus_to_data, dict_with_seq)
    assert df_out.loc[0, 'Location_in_sequence'] == 0
    assert df_out.loc[0, 'Transcript'] is None or isinstance(df_out.loc[0, 'Transcript'], str)

def test_get_id_from_annotations(tmp_path):
    content = """chr1\tHAVANA\tgene\t1000\t5000\t.\t+\t.\tgene_id \"GENE0001\"; gene_name \"GENE1\";
chr1\tHAVANA\tgene\t6000\t9000\t.\t+\t.\tgene_id \"GENE0002\"; gene_name \"GENE2\";
"""
    gtf_file = tmp_path / "test.gtf"
    gtf_file.write_text(content)

    from alternative_splicing import get_id_from_annotations
    gene_list = ["GENE1"]
    result = get_id_from_annotations(str(gtf_file), gene_list)
    assert result == {"GENE1": "GENE0001"}

def test_get_seq_for_genes():
    dict_with_id = {"GENE1": "ENSG00001"}
    locus_to_data = {"GENE1": DummyGene("ACGTACGT")}
    df = get_seq_for_genes(dict_with_id, locus_to_data)
    assert df.loc[0, "gene_sequence"] == "ACGTACGT"
