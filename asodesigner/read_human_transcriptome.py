from asodesigner.experiment import get_experiment
from asodesigner.file_utils import read_human_transcriptome_fasta_dict
from asodesigner.process_utils import run_off_target_hybridization_analysis


def main():
    fasta_dict = read_human_transcriptome_fasta_dict()

    experiment = get_experiment('SecondScrambled')

    print(experiment.target_sequence)

    run_off_target_hybridization_analysis(experiment, fasta_dict, organism='human')


if __name__ == "__main__":
    main()
