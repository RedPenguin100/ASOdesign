from asodesigner.experiment import get_experiment, maybe_create_experiment_folders, get_experiments
from asodesigner.file_utils import read_human_transcriptome_fasta_dict
from asodesigner.process_utils import run_off_target_hybridization_analysis, run_off_target_wc_analysis


def main():
    organism = 'human'

    fasta_dict = read_human_transcriptome_fasta_dict()
    # maybe_create_experiment_folders(this_experiment)
    experiments = get_experiments(['EntirePositiveControl'])

    for experiment in experiments:
        print(experiment.target_sequence)

        run_off_target_wc_analysis(experiment, fasta_dict, organism=organism)
        run_off_target_hybridization_analysis(experiment, fasta_dict, organism=organism)


if __name__ == "__main__":
    main()
