from Bio import SeqIO

from asodesigner.consts import HUMAN_TRANSCRIPTS_FASTA
from asodesigner.experiment import Experiment
from asodesigner.process_utils import run_off_target_hybridization_analysis
from asodesigner.target_finder import get_gfp_second_exp
from asodesigner.timer import Timer
from asodesigner.util import get_antisense


def main():
    with Timer() as t:
        with open(str(HUMAN_TRANSCRIPTS_FASTA), 'r') as handle:
            fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    print(f"Reading FASTA took: {t.elapsed_time}s")

    experiment = Experiment()
    experiment.name = "SecondScrambled"
    experiment.target_sequence = get_antisense(get_gfp_second_exp())
    experiment.l_values = [16, 17, 18, 19, 20, 21, 22]

    print(experiment.target_sequence)

    run_off_target_hybridization_analysis(experiment, fasta_dict)


if __name__ == "__main__":
    main()
