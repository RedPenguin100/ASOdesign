from asodesigner.consts import HUMAN_TRANSCRIPTS_FASTA, HUMAN_TRANSCRIPTS_FASTA_GZ
from asodesigner.experiment import get_experiments
from asodesigner.file_utils import get_fasta_dict_from_path
from asodesigner.process_utils import run_off_target_hybridization_analysis


def read_human_transcriptome_fasta_dict():
    if HUMAN_TRANSCRIPTS_FASTA.is_file():
        return get_fasta_dict_from_path(HUMAN_TRANSCRIPTS_FASTA)

    if HUMAN_TRANSCRIPTS_FASTA_GZ.is_file():
        return get_fasta_dict_from_path(HUMAN_TRANSCRIPTS_FASTA_GZ)

    raise FileNotFoundError(
        f"Did not find {HUMAN_TRANSCRIPTS_FASTA} or {HUMAN_TRANSCRIPTS_FASTA_GZ}, please consider the README.md")


def main():
    fasta_dict = read_human_transcriptome_fasta_dict()

    experiments = get_experiments(['SecondScrambled'])
    experiment = experiments[0]

    print(experiment.target_sequence)

    run_off_target_hybridization_analysis(experiment, fasta_dict, organism='human')


if __name__ == "__main__":
    main()
