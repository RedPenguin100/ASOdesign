from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import ndcg_score
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

# Spliting scheme
seed = 42
def splitting(df, split = 'A'):
    # Split randomly 20% of lines
    if split == 'A':
        train, test = train_test_split(df, test_size=0.2, random_state=seed)

    # Split by cell line: 20% of cell lines for test, 80% for train
    if split == 'B':
        cell_lines = df['Cell_line'].unique()
        cell_lines_train, cell_lines_test = train_test_split(cell_lines, test_size=0.2, random_state=seed)
        train = df[df['Cell_line'].isin(cell_lines_train)].copy()
        test = df[df['Cell_line'].isin(cell_lines_test)].copy()

    # Split 20% of each Cell line as test, rest as train
    if split == 'C':
        train_list = []
        test_list = []
        for cell_line, group in df.groupby('Cell_line'):
            group_train, group_test = train_test_split(group, test_size=0.2, random_state=seed)
            train_list.append(group_train)
            test_list.append(group_test)
        train = pd.concat(train_list).copy()
        test = pd.concat(test_list).copy()
        
    return train, test

def mard_per_cell_line(df, pred_col='predicted_score', true_col='Inhibition(%)',
                       group_col='Cell_line', k=5):
    """
    Compute MARD@all and MARD@k per cell line.
    """
    results = []

    for cell_line, g in df.groupby(group_col):
        g = g[[true_col, pred_col]].dropna().copy()
        n = len(g)
        if n < 2:
            continue

        k_eff = max(1, min(k, n))

        # ranks within this cell line only
        r_true = g[true_col].rank(ascending=False)
        r_pred = g[pred_col].rank(ascending=False)

        abs_diff = (r_pred - r_true).abs()

        # MARD over all items
        m_all = float(abs_diff.mean())

        # Top-k by truth/pred (deterministic tie-break)
        g["_ord"] = np.arange(n)
        top_true_idx = g.sort_values([true_col, "_ord"], ascending=[False, True]).head(k_eff).index
        top_pred_idx = g.sort_values([pred_col, "_ord"],  ascending=[False, True]).head(k_eff).index

        # MARD over items that are top-k by prediction
        diffs_pred = (r_pred.loc[top_pred_idx] - r_true.loc[top_pred_idx]).abs()
        m_k_pred = float(diffs_pred.mean())

        # Top-k hit rate (overlap)
        hit_rate = len(set(top_true_idx) & set(top_pred_idx)) / k_eff

        results.append({
            group_col: cell_line,
            'MARD@all': m_all,
            f'MARD@{k_eff}': m_k_pred,
            f'Top-{k_eff} hit rate': hit_rate,
        })

    return pd.DataFrame(results).sort_values(f'MARD@{k_eff}', ascending=False)

def ndcg_per_cell_line(df, pred_col='predicted_score', true_col='Inhibition(%)',
                       group_col='Cell_line', k=5):
    """
    Compute NDCG@all and NDCG@k per cell line.
    """
    results = []

    for cell_line, g in df.groupby(group_col):
        n_items = len(g)
        if n_items < 2:
            continue

        # Prepare arrays
        true_vals = g[true_col].to_numpy().reshape(1, -1)
        pred_vals = g[pred_col].to_numpy().reshape(1, -1)

        # Shift true values if any are negative (NDCG requires non-negative)
        if np.any(true_vals < 0):
            true_vals = true_vals - true_vals.min()

        # NDCG@all (no cutoff)
        ndcg_all = ndcg_score(true_vals, pred_vals)

        # NDCG@k (cutoff at k)
        ndcg_k = ndcg_score(true_vals, pred_vals, k=min(k, n_items))

        results.append({
            'Cell_line': cell_line,
            'Count': n_items,
            'NDCG@all': ndcg_all,
            f'NDCG@{k}': ndcg_k
        })
    return pd.DataFrame(results).sort_values(f'NDCG@{k}', ascending=False)

def permutation_corr_pvalue(x, y, method='pearson', n_permutations=10000, random_state=42):
    """Helper function: compute permutation-based correlation p-value."""
    rng = np.random.default_rng(random_state)

    if method == 'pearson':
        obs_corr, _ = pearsonr(x, y)
    elif method == 'spearman':
        obs_corr, _ = spearmanr(x, y)
    else:
        raise ValueError("method must be 'pearson' or 'spearman'")

    perm_corrs = []
    for _ in range(n_permutations):
        y_perm = rng.permutation(y)
        if method == 'pearson':
            r, _ = pearsonr(x, y_perm)
        else:
            r, _ = spearmanr(x, y_perm)
        perm_corrs.append(r)

    perm_corrs = np.array(perm_corrs)
    p_value_perm = np.mean(np.abs(perm_corrs) >= np.abs(obs_corr))

    return obs_corr, p_value_perm

def correlation_per_cell_line_with_permutation(df, pred_col='predicted_score', true_col='Inhibition(%)',
                                               group_col='Cell_line', n_permutations=10000, random_state=42):
    """
    Per-cell-line correlation metrics with both scipy p-values and permutation p-values.
    """
    results = []

    for cell_line, g in df.groupby(group_col):
        if len(g) < 3:
            continue  # need at least 3 samples for Pearson/Spearman

        y_true = g[true_col].to_numpy()
        y_pred = g[pred_col].to_numpy()

        # Pearson (SciPy)
        pear_corr, pear_p_scipy = pearsonr(y_true, y_pred)
        # Pearson (Permutation)
        pear_corr_perm, pear_p_perm = permutation_corr_pvalue(y_true, y_pred, method='pearson',
                                                              n_permutations=n_permutations,
                                                              random_state=random_state)

        # Spearman (SciPy)
        spear_corr, spear_p_scipy = spearmanr(y_true, y_pred)
        # Spearman (Permutation)
        spear_corr_perm, spear_p_perm = permutation_corr_pvalue(y_true, y_pred, method='spearman',
                                                                n_permutations=n_permutations,
                                                                random_state=random_state)

        results.append({
            'Cell_line': cell_line,
            'Count': len(g),
            'Pearson': pear_corr,
            'Pearson_p_scipy': pear_p_scipy,
            'Pearson_p_perm': pear_p_perm,
            'Spearman': spear_corr,
            'Spearman_p_scipy': spear_p_scipy,
            'Spearman_p_perm': spear_p_perm
        })

    return pd.DataFrame(results).sort_values('Pearson', ascending=False)
