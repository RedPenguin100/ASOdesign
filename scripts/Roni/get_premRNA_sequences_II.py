import os
import pandas as pd
from off_target_functions import dna_to_rna_reverse_complement, normalize_chrom
import pickle

# Get the directory of the current script
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# Full path to cell line data folder
DATA_DIR = os.path.join(PROJECT_ROOT, "scripts", "data_genertion", "cell_line_expression")

cell_line2dataII = {
    "ACH-001328": [os.path.join(DATA_DIR, "ACH-001328_transcriptome_premRNA.csv"), os.path.join(DATA_DIR, "ACH-001328_mutations.csv")],
    "ACH-000463": [os.path.join(DATA_DIR, "ACH-000463_transcriptome_premRNA.csv"), os.path.join(DATA_DIR, "ACH-000463_mutations.csv")],
    "ACH-001188": [os.path.join(DATA_DIR, "ACH-001188_transcriptome_premRNA.csv"), os.path.join(DATA_DIR, "ACH-001188_mutations.csv")],
    "ACH-001086": [os.path.join(DATA_DIR, "ACH-001086_transcriptome_premRNA.csv"), os.path.join(DATA_DIR, "ACH-001086_mutations.csv")],
    "ACH-000739": [os.path.join(DATA_DIR, "ACH-000739_transcriptome_premRNA.csv"), os.path.join(DATA_DIR, "ACH-000739_mutations.csv")],
    "ACH-000232": [os.path.join(DATA_DIR, "ACH-000232_transcriptome_premRNA.csv"), os.path.join(DATA_DIR, "ACH-000232_mutations.csv")]
}

def all_premRNA_sequences(cell_line_dict, genome_pkl, annotation_pkl):
    with open(annotation_pkl, 'rb') as f:
        gtf_dict = pickle.load(f)
    print('Read GTF.')
    with open(genome_pkl, 'rb') as f:
        fasta_dict = pickle.load(f)
    print('Read FASTA.')
    print(f'fasta keys: {fasta_dict.keys()}')

    all_results = {}
    for cell_line, [exp_data, mut_data] in cell_line_dict.items():
        print(f"Processing {cell_line}...")

        exp_data = pd.read_csv(exp_data)
        for idx, row in exp_data.iterrows():

            if not pd.isna(row["Original Transcript Sequence"]):
                continue
            else:
                transcript_id_full = row["Transcript_ID"]
                transcript_id = transcript_id_full.split('.')[0]  # remove version

                transcript_info = gtf_dict.get(transcript_id_full) or gtf_dict.get(transcript_id)

                if not transcript_info:
                    continue

                chrom = transcript_info["chrom"]
                strand = transcript_info["strand"]
                start = transcript_info["start"]
                end = transcript_info["end"]

                try:
                    seq = fasta_dict[chrom][start - 1:end]  # extract using 0-based indexing
                    if strand == "-":
                        seq = dna_to_rna_reverse_complement(seq)
                    else:
                        seq = seq.replace("T", "U")
                    exp_data.at[idx, 'Original Transcript Sequence'] = seq
                except Exception as e:
                    print(chrom in fasta_dict)  # Should return True or False
                    print("Trying chrom:", chrom)
                    print(f"Error extracting for {transcript_id_full}: {e}")

        output_path = f"{cell_line}_transcriptome_premRNA.csv"
        exp_data.to_csv(output_path, index=False)
        print(f"Saved {output_path}")

        all_results[cell_line] = exp_data

    return all_results


genome_path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/_whole_genomic_sequence.pkl'
annotation_path = '/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/_gtf_annotations.pkl'
all_premRNA_sequences(cell_line2dataII, genome_path, annotation_path)