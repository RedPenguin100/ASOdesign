import unittest
import pandas as pd
from extracting_data_from_article import (
    find_if_seq_in_gene,
    load_csv,
    filter_by_mod,
    read_fasta_biopython,
    filter_cell_line_human
)


class TestFindSeqInGene(unittest.TestCase):
    def test_sequence_found(self):
        self.assertEqual(find_if_seq_in_gene("ATG", "CCATGGA"), 2)

    def test_sequence_not_found(self):
        self.assertEqual(find_if_seq_in_gene("TGA", "CCATGGA"), -1)

    def test_case_insensitivity(self):
        self.assertEqual(find_if_seq_in_gene("atg", "CCATGGA"), 2)

    def test_full_match(self):
        self.assertEqual(find_if_seq_in_gene("CCATGGA", "CCATGGA"), 0)

    def test_empty_sequence(self):
        self.assertEqual(find_if_seq_in_gene("", "CCATGGA"), 0)

    def test_empty_gene(self):
        self.assertEqual(find_if_seq_in_gene("ATG", ""), -1)

    def test_longer_sequence_than_gene(self):
        self.assertEqual(find_if_seq_in_gene("CCATGGAAT", "CCATGGA"), -1)


class TestFilterByMod(unittest.TestCase):
    def setUp(self):
        self.test_df = pd.DataFrame({
            "Modification": ["Mod1", "Mod2", "Mod1", "Mod3"],
            "Inhibition(%)": [50, 30, 80, 60]
        })

    def test_filter_modification(self):
        result = filter_by_mod(self.test_df, "Mod1")
        self.assertTrue(all(result["Modification"] == "Mod1"))
        self.assertEqual(result["Inhibition(%)"].tolist(), [80, 50])  # Sorted


class TestFilterCellLineHuman(unittest.TestCase):
    def setUp(self):
        # Test data for the main dataframe
        self.test_df = pd.DataFrame({"Cell_line": ["A", "B", "C"]})

        # Test data for the cell line info dataframe
        self.cell_line_info = pd.DataFrame({
            "Cell_line": ["A", "B", "C"],
            "Cell line organism": ["human", "mouse", "human"]
        })

    def test_filter_human_cell_lines(self):
        # Test the filter function with the provided cell_line_info dataframe
        filtered_df = filter_cell_line_human(self.test_df, self.cell_line_info)

        # Check if only human cell lines are returned
        self.assertSetEqual(set(filtered_df["Cell_line"]), {"A", "C"})


if __name__ == "__main__":
    unittest.main()
