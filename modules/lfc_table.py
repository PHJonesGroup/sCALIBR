import pandas as pd
import numpy as np

def lfc_table(tab, baseline, cond1):
    """
   Build a per-replicate log2 fold-change table for a set of gRNAs.

    Parameters
    ----------
    tab : pandas.DataFrame
        DataFrame [gRNA, gene, <normalised count columns>]
    baseline : str
        Baseline condition name / column prefix
    cond1 : str
        Treatment condition name / column prefix

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

    if len(baseline) != len(cond1):
        raise ValueError(f"baseline ({baseline}) and cond1 ({cond1}) must be the same length")

    suffixes = [_pairing_suffix(b, c) for b, c in zip(baseline, cond1)]
    dup = sorted({s for s in suffixes if suffixes.count(s) > 1})
    if dup:
        raise ValueError(f"Replicate labels {dup} are not unique for baseline {baseline}; "
                         f"two replicates would overwrite each other's LFC column")

    for b, c, suffix in zip(baseline, cond1, suffixes):
        if b not in tab.columns or c not in tab.columns:
            raise ValueError(f"Missing '{b}' or '{c}'. Columns: {tab.columns.tolist()}")
        out[f"lfc_{suffix}"] = np.log2(
            (tab[c].astype(float) + eps) / (tab[b].astype(float) + eps)
        ).values

    return out

def _pairing_suffix(baseline_col, cond1_col):
    """Derive a short label for a baseline/cond1 pairing, e.g. ('T0_F','T1_F+') -> 'F',
    ('T0_F_R1','T1_F_R1') -> 'F_R1'."""
    # take everything after the first underscore of the baseline column (drops the
    # time-point prefix), stripped of +/-
    tail = baseline_col.split('_', 1)[-1]           # 'T0_F_R1' -> 'F_R1'
    return tail.strip('+-') or baseline_col