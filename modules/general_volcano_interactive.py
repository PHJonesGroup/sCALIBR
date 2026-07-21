import numpy as np
import pandas as pd

def general_volcano_interactive(alf, sfdr_corr, thr_scoreh, thr_scored, score, fdr, cond, genes):
    thr_fdr = alf
    fdr_corr = fdr + sfdr_corr  # avoid log10(0)
    LPV = -np.log10(fdr_corr)

    T_gene = pd.DataFrame({'gene': genes, 'score': score, 'fdr_corr': fdr_corr, 'LPV': LPV})

    indh = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s >= thr_scoreh)]
    indd = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s <= thr_scored)]

    T_hs = T_gene.loc[indh].copy() if indh else pd.DataFrame()
    T_ds = T_gene.loc[indd].copy() if indd else pd.DataFrame()

    return LPV, indh, indd, T_ds, T_hs, T_gene
