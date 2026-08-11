import pandas as pd

def norm_table_individ(norm_dat, baseline, cond1):
    """
    Split normalized data into one table per condition (baseline, cond1),
    keeping the sgRNA and gene ID columns plus that condition's replicate columns.
    
    Parameters
    ----------
    norm_dat : pandas.DataFrame
        Normalised counts; ID columns first, then '<condition>_<rep>' count columns.
    baseline : str
        Baseline condition name / column prefix
    cond1: str
        Treatment condition name / column prefix

    Returns
    -------
    tab_norm_baseline : pandas.DataFrame
        [sgRNA_name, gene, <baseline replicate columns>].
    tab_norm_cond1 : pandas.DataFrame
        [sgRNA_name, gene, <cond1 replicate columns>].
    
    Raises
    ------
    ValueError
        If a required count column (baseline/cond1) is missing from ``norm_dat``
    """
    def build_table(cond):
        # all data columns belonging to this condition
        cond_cols = [c for c in norm_dat.columns if c.startswith(f"{cond}")]
        if not cond_cols:
            raise ValueError(f"No columns found for condition '{cond}'")

        tab = pd.DataFrame({
            'sgRNA_name': norm_dat.iloc[:, 0].values,
            'gene':       norm_dat.iloc[:, 1].values,
            'gene_type':  norm_dat.iloc[:, 2].values,
        })
        for c in cond_cols:
            tab[c] = norm_dat[c].values
        return tab

    tab_norm_baseline = build_table(baseline)
    tab_norm_cond1 = build_table(cond1)

    return tab_norm_baseline, tab_norm_cond1