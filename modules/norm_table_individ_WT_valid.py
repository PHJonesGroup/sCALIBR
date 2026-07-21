import pandas as pd

def norm_table_individ_WT_valid(norm_dat, raw_ind, cond1, cond2, rep=None):
    """
    Split normalized data into one table per condition (cond1, cond2),
    keeping the sgRNA and gene ID columns plus that condition's replicate columns.

    Inputs:
        norm_dat: pandas DataFrame — ID columns first, then data columns
        raw_ind:  pandas DataFrame with sgRNA names (col 0) and gene names (col 1)
        cond1, cond2: condition name prefixes

    Outputs:
        tab_norm_cond1, tab_norm_cond2: DataFrames with
            ['sgRNA_name', 'gene', <replicate columns for that condition>]
    """
    id_cols = ['sgRNA_name', 'gene']

    def build_table(cond):
        # all data columns belonging to this condition, in original order
        cond_cols = [c for c in norm_dat.columns if c.startswith(f"{cond}_")]
        if not cond_cols:
            raise ValueError(f"No columns found for condition '{cond}_'")

        tab = pd.DataFrame({
            'sgRNA_name': raw_ind.iloc[:, 0].values,
            'gene':       raw_ind.iloc[:, 1].values,
        })
        # attach each replicate column under its own name
        for c in cond_cols:
            tab[c] = norm_dat[c].values
        return tab

    tab_norm_cond1 = build_table(cond1)
    tab_norm_cond2 = build_table(cond2)

    return tab_norm_cond1, tab_norm_cond2