import pandas as pd
import numpy as np

def lfc_table(tab, baseline, cond1, rep_pairs):
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
    rep_pairs : list of str or None
        Name of replicates

    Returns
    -------
    out : pandas.DataFrame
        [gRNA, gene, lfc]                 if no replicates, or
        [gRNA, gene, lfc_<rep>, ...]      one LFC column per replicate.

    Raises
    ------
    ValueError
        If a required count column (baseline/cond1, with rep suffix if applicable)
        is missing from ``tab``
    """
    eps = 1e-6
    out = pd.DataFrame({
        'gRNA': tab.iloc[:, 0].values,
        'gene': tab.iloc[:, 1].astype(str).values,
        'gene_type' : tab.iloc[:, 2].astype(str).values,
    })

    if not rep_pairs:                      # None or empty -> no replicate suffixes
        c1, c2 = baseline, cond1
        if c1 not in tab.columns or c2 not in tab.columns:
            raise ValueError(f"Missing {c1} or {c2}. Columns: {tab.columns.tolist()}")
        out['lfc'] = np.log2(
            (tab[c2].astype(float) + eps) / (tab[c1].astype(float) + eps)
        ).values
    else:                                  # replicate suffixes present
        for r in rep_pairs:
            c1, c2 = f"{baseline}_{r}", f"{cond1}_{r}"
            if c1 not in tab.columns or c2 not in tab.columns:
                raise ValueError(f"Missing {c1} or {c2}")
            out[f"lfc_{r}"] = np.log2(
                (tab[c2].astype(float) + eps) / (tab[c1].astype(float) + eps)
            ).values
    return out
