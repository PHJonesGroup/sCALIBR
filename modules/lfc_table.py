import pandas as pd
import numpy as np

def lfc_table(tab, baseline, cond1):
    """
   Build a per-replicate log2 fold-change table for a set of gRNAs.

    Parameters
    ----------
    tab : pandas.DataFrame
        DataFrame [gRNA, gene, <normalised count columns>]
    baseline : str or list of str
        Baseline sample column(s), paired in order with cond1; a single
        baseline is shared by every cond1 replicate.
    cond1 : str or list of str
        Treatment sample column(s).

    Returns
    -------
    out : pandas.DataFrame
        [gRNA, gene, lfc]                 if no replicates, or
        [gRNA, gene, lfc_<rep>, ...]      one LFC column per replicate.

    Raises
    ------
    ValueError
        If a required count column (baseline/cond1, with rep suffix if applicable) is missing from tab
    """
    eps = 1e-6
    out = pd.DataFrame({
        'gRNA':      tab.iloc[:, 0].values,
        'gene':      tab.iloc[:, 1].astype(str).values,
        'gene_type': tab.iloc[:, 2].astype(str).values,
    })

    # normalise to lists
    if isinstance(baseline, str):
        baseline = [baseline]
    if isinstance(cond1, str):
        cond1 = [cond1]

    # one baseline shared by every treatment replicate (e.g. a single T0 or plasmid)
    if len(baseline) == 1 and len(cond1) > 1:
        baseline = baseline * len(cond1)
    if len(baseline) != len(cond1):
        raise ValueError(f"baseline ({baseline}) and cond1 ({cond1}) must be the same length, "
                         f"or baseline a single sample")

    # label replicates by baseline, or by treatment when baselines are shared
    shared = len(set(baseline)) < len(baseline)
    suffixes = _replicate_labels(cond1 if shared else baseline)

    for b, c, suffix in zip(baseline, cond1, suffixes):
        if b not in tab.columns or c not in tab.columns:
            raise ValueError(f"Missing '{b}' or '{c}'. Columns: {tab.columns.tolist()}")
        out[f"lfc_{suffix}"] = np.log2(
            (tab[c].astype(float) + eps) / (tab[b].astype(float) + eps)
        ).values

    return out

def _replicate_labels(baseline):
    """One unique label per sample column, used as the LFC column suffix.

    Tries, in order, the first rule that gives unique labels:
    1. _pairing_suffix: drop the time-point prefix ('T0_F_R1' -> 'F_R1');
    2. drop the '_' parts shared by every baseline column
       (['F_T0_rep1', 'M_T0_rep1'] -> ['F_rep1', 'M_rep1']);
    3. the full baseline column name.
    """
    def unique(labels):
        return all(labels) and len(set(labels)) == len(labels)

    labels = [_pairing_suffix(b) for b in baseline]
    if unique(labels):
        return labels

    parts = [b.split('_') for b in baseline]
    shared = set.intersection(*(set(p) for p in parts))
    labels = ['_'.join(t for t in p if t not in shared).strip('+-') for p in parts]
    if unique(labels):
        return labels

    labels = list(baseline)
    if not unique(labels):
        raise ValueError(f"Baseline columns {baseline} are not unique; "
                         f"two replicates would overwrite each other's LFC column")
    return labels

def _pairing_suffix(baseline_col):
    """Derive a short replicate label from the baseline column, e.g. 'T0_F' -> 'F',
    'T0_F_R1' -> 'F_R1'."""
    # take everything after the first underscore of the baseline column (drops the
    # time-point prefix), stripped of +/-
    tail = baseline_col.split('_', 1)[-1]           # 'T0_F_R1' -> 'F_R1'
    return tail.strip('+-') or baseline_col