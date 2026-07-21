import numpy as np
import pandas as pd

import numpy as np
import pandas as pd

def filter_pattern_distri(raw_ind, pat, n, st, en, step, cond1, cond2, rep_pairs):
    # 1. Extract genes and gRNAs
    genes = raw_ind.iloc[:, 1].astype(str)
    gRNA  = raw_ind.iloc[:, 0]

    # 2. Pattern filter on first n characters
    matches = genes.apply(lambda g: g[:n].lower() == pat.lower())
    ind_out = matches[matches].index
    ind_in  = matches[~matches].index

    indiv_Chr   = raw_ind.loc[ind_out]
    indiv_noChr = raw_ind.loc[ind_in]
    num_out_in  = len(ind_out)

    # 3. One LFC per replicate pair: log2(cond1_<rep> / cond2_<rep>)
    eps = 1e-6
    LFC_cols = {}
    for r in rep_pairs:
        c1 = f"{cond1}_{r}"
        c2 = f"{cond2}_{r}"
        if c1 not in raw_ind.columns or c2 not in raw_ind.columns:
            raise ValueError(f"Missing {c1} or {c2}")
        if len(indiv_Chr):
            LFC_cols[f"lfc_{r}"] = np.log2(
                (indiv_Chr[c1].astype(float) + eps) /
                (indiv_Chr[c2].astype(float) + eps)
            ).values
        else:
            LFC_cols[f"lfc_{r}"] = np.empty(0)

    LFC = np.column_stack(list(LFC_cols.values())) if LFC_cols else np.empty((len(indiv_Chr), 0))

    # 4. Histogram of pooled LFC across all replicate pairs (for the FM summary)
    bins = np.arange(st, en + step, step)
    flat = LFC[np.isfinite(LFC)]
    hissFM, _ = np.histogram(flat, bins)
    percFM = 100 * hissFM / hissFM.sum() if hissFM.sum() > 0 else np.zeros_like(hissFM, dtype=float)

    all_cols = [f"{cond1}_{r}" for r in rep_pairs] + [f"{cond2}_{r}" for r in rep_pairs]
    s_chr = indiv_Chr[all_cols].astype(float).sum().values
    st1, en1 = st, en

    # 5. LFC table: ids + one lfc column per replicate pair
    T_lfc_pat = pd.DataFrame({
        'gRNA':  gRNA.loc[ind_out].values,
        'gene':  genes.loc[ind_out].values,     # 'gene' (singular) for consistency
    })
    for name, vals in LFC_cols.items():
        T_lfc_pat[name] = vals

    return (
        indiv_Chr, indiv_noChr,
        ind_out.tolist(), ind_in.tolist(), num_out_in,
        bins, hissFM, percFM,
        LFC,
        s_chr, st1, en1, T_lfc_pat
    )