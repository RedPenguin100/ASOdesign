import RNA
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from access_calculator import AccessCalculator

def calculate_average_mfe_per_nucleotide(sequence, window_size=120, steps=[6,12,18]):

    sequence_length = len(sequence)
    energy_values = np.zeros((len(steps), sequence_length))
    counts = np.zeros((len(steps), sequence_length))

    for step_idx, step in enumerate(steps):
        for i in range(0, sequence_length - window_size + 1, step):
            subseq = sequence[i:i + window_size]
            _, mfe = RNA.fold(subseq)
            mfe_per_nt = -mfe / window_size

            # Assign energy value to all positions in the current window
            for j in range(i, i + window_size):
                energy_values[step_idx, j] += mfe_per_nt
                counts[step_idx, j] += 1

    counts[counts == 0] = 1 # Avoid division by zero
    avg_energies_per_step = energy_values / counts # Element-wise division
    avg_energies_nt = np.mean(avg_energies_per_step, axis=0)

    return avg_energies_nt

def sliding_window_sum(lst, k):
    """"
    example: ([1,2,3,4,5,6],3) -> [6,9,12,15]
    """
    return np.array([sum(lst[i:i+k]) for i in range(len(lst) - k + 1)])

def show_normalized_plot(energies, title):
    x_min, x_max = energies.min(), energies.max()
    x_normalized = (energies - x_min) / (x_max - x_min)

    plt.plot(x_normalized, marker='o', linestyle='-', color='b')
    plt.xlabel('position')
    plt.ylabel('normalized value')
    plt.title(title)
    plt.grid(True)
    plt.show()

def spearman_corr(s1, s2):
    min_len = min(len(s1), len(s2))
    return sp.stats.spearmanr(s1[:min_len], s2[:min_len])

if __name__ == "__main__":
    YEAST_M_CHERRY = 'ATGTCTAAGGGGGAAGAAGACAATATGGCGATTATTAAAGAGTTTATGAGATTTAAAGTACATATGGAAGGAAGTGTTAATGGTCACGAGTTTGAGATCGAAGGTGAAGGTGAAGGTCGTCCATATGAGGGTACGCAAACAGCAAAACTAAAGGTGACTAAAGGGGGACCATTACCTTTCGCTTGGGATATACTGTCACCACAATTCATGTACGGATCGAAAGCTTACGTAAAGCACCCGGCCGACATTCCTGATTATTTAAAGTTGTCTTTCCCTGAAGGGTTCAAATGGGAAAGAGTTATGAATTTTGAGGATGGAGGTGTTGTGACGGTAACTCAAGATTCATCTTTGCAAGATGGCGAATTCATTTATAAAGTTAAATTGAGAGGAACTAACTTTCCAAGCGATGGTCCAGTCATGCAAAAAAAGACCATGGGCTGGGAAGCTAGCTCAGAACGGATGTACCCGGAAGACGGCGCATTAAAGGGAGAGATCAAGCAGCGACTTAAGTTAAAAGATGGCGGGCATTATGATGCAGAAGTAAAGACAACCTACAAAGCCAAAAAACCCGTGCAGCTGCCTGGTGCGTATAATGTTAACATAAAACTAGACATTACATCCCACAACGAAGACTACACTATAGTCGAACAATACGAAAGGGCAGAAGGTAGACATTCGACAGGTGGTATGGATGAGTTGTATAAATAA'.replace('T', 'U')
    RNA_sequence = YEAST_M_CHERRY

    avg_energies = calculate_average_mfe_per_nucleotide(RNA_sequence, window_size=45, steps=[6,12,18])

    show_normalized_plot(avg_energies, 'Average mfe per nt')
    print(len(avg_energies))
    avg_energies_6 = sliding_window_sum(avg_energies, 6)
    show_normalized_plot(avg_energies_6, 'Average MFE for 6 nucleotides segments')
    print(len(avg_energies_6))

    g_seq = YEAST_M_CHERRY
    g_access_size = 12
    g_min_gc = 0
    g_max_gc = 100
    g_gc_ranges = 1
    # g_temperature = 37
    g_access_win_size = 120
    g_access_seed_size = 3
    g_access_seed_sizes = [g_access_seed_size * m for m in range(1, 4)]
    g_df = AccessCalculator.calc(
        g_seq, g_access_size, g_min_gc, g_max_gc, g_gc_ranges, g_access_win_size, g_access_seed_sizes)

    energies_yehuda_3_avg = g_df['3_avg']
    energies_yehuda_6_avg = g_df['6_avg']
    energies_yehuda_avg_access = g_df['avg_access']

    show_normalized_plot(energies_yehuda_3_avg, 'energies_yehuda_3_avg')
    show_normalized_plot(energies_yehuda_6_avg, 'energies_yehuda_6_avg')
    show_normalized_plot(energies_yehuda_avg_access, 'energies_yehuda_avg_access')

    print(spearman_corr(avg_energies_6, energies_yehuda_6_avg))
    print(spearman_corr(avg_energies_6, energies_yehuda_avg_access))



