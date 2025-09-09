import numpy as np
import ViennaRNA


def get_sense_with_flanks(pre_mrna: str, sense_start: int, sense_length: int, flank_size: int) -> str:
    """
    Re  turns the sense sequence with `flank_size` nucleotides on each side (if available).
    If near the edge, it will not go out of bounds.

    Parameters:
    - pre_mrna: The full pre-mRNA sequence (5' -> 3')
    - sense_start: Start index of the sense sequence within pre_mrna
    - sense_length: Length of the sense sequence (usually same as antisense length)
    - flank_size: Number of nucleotides to include on each side (upstream and downstream)

    Returns:
    - str: The flanked sense sequence
    """
    # Ensure indices are within bounds
    start = max(0, sense_start - flank_size)
    end = min(len(pre_mrna), sense_start + sense_length + flank_size)

    return pre_mrna[start:end]

def calculate_avg_mfe_over_sense_region(sequence, sense_start, sense_length, flank_size=120, window_size=45, step=7):
    sequence = str(sequence).upper().replace('T', 'U')
    sequence_length = len(sequence)
    energy_values = np.zeros(sequence_length)
    counts = np.zeros(sequence_length)

    for i in range(0, sequence_length - window_size + 1, step):
        subseq = sequence[i:i + window_size]
        fc = ViennaRNA.fold_compound(subseq)
        _, mfe = fc.mfe()
        mfe_per_nt = mfe / window_size

        for j in range(i, i + window_size):
            energy_values[j] += mfe_per_nt
            counts[j] += 1

    counts[counts == 0] = 1
    avg_energies = energy_values / counts

    flank_start = max(0, sense_start - flank_size)
    sense_start_in_flank = sense_start - flank_start
    sense_end_in_flank = sense_start_in_flank + sense_length

    if 0 <= sense_start_in_flank < sequence_length and sense_end_in_flank <= sequence_length:
        return np.mean(avg_energies[sense_start_in_flank:sense_end_in_flank])
    else:
        return np.nan

def calculate_mfe_over_edges_sense_region(sequence, sense_start, sense_length, flank_size=45, window_size=45, step=7):
    sequence = str(sequence).upper().replace('T', 'U')
    sequence_length = len(sequence)
    energy_values = np.zeros(sequence_length)
    counts = np.zeros(sequence_length)

    for i in range(0, sequence_length - window_size + 1, step):
        subseq = sequence[i:i + window_size]
        fc = ViennaRNA.fold_compound(subseq)
        _, mfe = fc.mfe()
        mfe_per_nt = mfe / window_size

        for j in range(i, i + window_size):
            energy_values[j] += mfe_per_nt
            counts[j] += 1

    counts[counts == 0] = 1
    avg_energies = energy_values / counts

    flank_start = max(0, sense_start - flank_size)
    sense_start_in_flank = sense_start - flank_start
    sense_end_in_flank = sense_start_in_flank + sense_length

    if 0 <= sense_start_in_flank < sequence_length and sense_end_in_flank <= sequence_length:
        return np.mean(np.concatenate([(avg_energies[sense_start_in_flank:sense_start_in_flank+4]), (avg_energies[sense_end_in_flank-4:sense_end_in_flank])]))
    else:
        return np.nan

def calculate_min_mfe_over_sense_region(sequence, sense_start, sense_length, flank_size=120, window_size=45, step=7):
    sequence = str(sequence).upper().replace('T', 'U')
    sequence_length = len(sequence)
    energy_values = np.zeros(sequence_length)
    counts = np.zeros(sequence_length)

    for i in range(0, sequence_length - window_size + 1, step):
        subseq = sequence[i:i + window_size]
        fc = ViennaRNA.fold_compound(subseq)
        _, mfe = fc.mfe()
        mfe_per_nt = mfe / window_size

        for j in range(i, i + window_size):
            energy_values[j] += mfe_per_nt
            counts[j] += 1

    counts[counts == 0] = 1
    avg_energies = energy_values / counts

    flank_start = max(0, sense_start - flank_size)
    sense_start_in_flank = sense_start - flank_start
    sense_end_in_flank = sense_start_in_flank + sense_length

    if 0 <= sense_start_in_flank < sequence_length and sense_end_in_flank <= sequence_length:
        return np.min(avg_energies[sense_start_in_flank:sense_end_in_flank])
    else:
        return np.nan

