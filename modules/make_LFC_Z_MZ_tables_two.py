import pandas as pd
import numpy as np

def make_LFC_Z_MZ_tables_two(zGE, MZ, Z):
    lfc_cols = list(zGE.columns[2:])
    n_reps = len(lfc_cols)
    
    return pd.DataFrame({
        "gRNA": np.repeat(zGE.iloc[:, 0].values, n_reps),
        "gene": np.repeat(zGE.iloc[:, 1].values, n_reps),
        "rep":  np.tile(lfc_cols, len(zGE)),
        "LFC":  zGE[lfc_cols].astype(float).values.ravel(),
        "Z":    np.asarray(Z, dtype=float),
        "MZ":   np.asarray(MZ, dtype=float),
    })

