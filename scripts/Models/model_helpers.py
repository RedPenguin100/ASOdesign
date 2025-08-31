import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scripts.data_genertion.consts import *


def correction(df):
    return df[VOLUME] / (df[VOLUME] + 10)

def log_inhibition_to_regular(log_inhibition, log_correction):
    return 100 * (-np.exp(-log_inhibition) + log_correction)

def get_threshold_for_percentile(percentile, df):
    k = max(1, int(np.ceil((percentile / 100) * len(df))))  # desired count
    threshold = np.partition(df, -k)[-k]  # value at 99th percentile
    return threshold


def evaluate_top(model, test, metric, features, log_correction, top_k=None, plot=False):
    test_filtered = test.copy()
    test_pred = model.predict(test_filtered[features].to_numpy())

    if top_k is not None:
        test_mask = np.zeros_like(test_pred, dtype=bool)
        top_k_idx = np.argsort(test_pred)[-top_k:]
        test_mask[top_k_idx] = True
    else:
        percentile = 0.5
        test_threshold = get_threshold_for_percentile(percentile, test_pred)
        test_mask = test_pred > test_threshold
        top_k = test_mask.sum()

    y_test = test_filtered[metric].to_numpy()
    if metric == INHIBITION:
        top_test = y_test[test_mask]
        all_inhib = pd.Series(y_test)
    elif metric == 'log_inhibition':
        top_test = log_inhibition_to_regular(y_test[test_mask], log_correction)
        all_inhib = pd.Series(log_inhibition_to_regular(y_test, log_correction))
    else:
        top_test = log_inhibition_to_regular(y_test[test_mask] * correction(test_filtered[test_mask]), log_correction)
        all_inhib = log_inhibition_to_regular(y_test * correction(test_filtered), log_correction)


    # ---------- best possible ----------
    top_best = all_inhib.nlargest(top_k).to_numpy()  # no index confusion

    # random sample via mask
    # k already computed as the desired sample size
    rand_mask = np.zeros(len(y_test), dtype=bool)
    rand_mask[np.random.choice(len(y_test), top_k, replace=False)] = True

    top_rand_test = all_inhib[rand_mask]

    mean_model, mean_best, mean_random = np.mean(top_test), np.mean(top_best), np.mean(top_rand_test)

    if plot:
        bins = np.array([-50, 0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 97, 100])
        x = (bins[:-1] + bins[1:]) / 2

        fig, ax = plt.subplots(figsize=(4, 4))

        # --- Predicted: thick line, darker fill ---------------------------------
        h_pred, _ = np.histogram(top_test, bins=bins)
        ax.plot(x, h_pred, drawstyle='steps-mid', lw=3.2, color='C0',
                label='Predicted', zorder=5)
        ax.fill_between(x, 0, h_pred, step='mid', alpha=0.35, color='C0', zorder=4)

        # --- Others: thinner, lighter -------------------------------------------
        for arr, lab, c in [(top_rand_test, 'Random', 'C1'),
                            (top_best, 'Best', 'C2')]:
            h, _ = np.histogram(arr, bins=bins)
            ax.plot(x, h, drawstyle='steps-mid', lw=1.3, color=c, alpha=0.7, label=lab, zorder=3)
            ax.fill_between(x, 0, h, step='mid', alpha=0.10, color=c, zorder=2)

        # ax.set_yscale('log')        # drop if you don't want log
        ax.set_xlabel('Actual inhibition')
        ax.legend(frameon=False)
        fig.tight_layout()
        plt.show()

    return mean_model, mean_best, mean_random