import gzip

from pathlib import Path
from Bio import SeqIO
import warnings

from asodesigner.timer import Timer


def get_fasta_dict_from_path(fasta_path: Path):
    with Timer() as timer:
        if fasta_path.suffix == ".gz":
            warnings.warn(f"Fasta is compressed, consider decompressing for performance. To unzip, run gunzip {fasta_path}")
            with gzip.open(str(fasta_path), 'rt') as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        else:
            with open(str(fasta_path), 'r') as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    print(f"Time took to read fasta: {timer.elapsed_time}")
    return fasta_dict