import numpy as np
import pandas as pd 

def zGE_target_distri(genes_nnt, zGE_genes, T_norm, st, en, step):
    """
    Separate zGE and non-zGE target genes and compute LFC distributions.
    """
    zGE_set = set(zGE_genes)
    nzGE_set = set(T_norm.iloc[:, 1]) - zGE_set

    # Select rows for zGE genes
    T_zGE = T_norm[T_norm.iloc[:, 1].isin(zGE_set)].copy()
    # Select rows for non-zGE genes
    T_nzGE = T_norm[T_norm.iloc[:, 1].isin(nzGE_set)].copy()

    # Compute LFCs for zGE genes
    counts_zGE = T_zGE.iloc[:, 2:6].astype(float)
    LFC_1 = np.log2((counts_zGE['T1_F'] + 1e-6) / (counts_zGE['T0_F'] + 1e-6))
    LFC_2 = np.log2((counts_zGE['T1_M'] + 1e-6) / (counts_zGE['T0_M'] + 1e-6))
    LFC_sum = LFC_1 + LFC_2

    bins = np.arange(st, en + step, step)
    hissFM, _ = np.histogram(LFC_sum, bins)
    percFM = 100 * hissFM / hissFM.sum() if hissFM.sum() > 0 else np.zeros_like(hissFM)

    hiss1, _ = np.histogram(LFC_1, bins)
    hiss2, _ = np.histogram(LFC_2, bins)
    perc1 = 100 * hiss1 / hiss1.sum() if hiss1.sum() > 0 else np.zeros_like(hiss1)
    perc2 = 100 * hiss2 / hiss2.sum() if hiss2.sum() > 0 else np.zeros_like(hiss2)

    T_lfc_zGE = pd.DataFrame({
        'gRNA': T_zGE.iloc[:, 0].values,
        'gene': T_zGE.iloc[:, 1].values,
        'lfc1': LFC_1.values,
        'lfc2': LFC_2.values
    })

    # Compute LFCs for non-zGE genes
    counts_nzGE = T_nzGE.iloc[:, 2:6].astype(float)
    LFC_1t = np.log2((counts_nzGE['T1_F'] + 1e-6) / (counts_nzGE['T0_F'] + 1e-6))
    LFC_2t = np.log2((counts_nzGE['T1_M'] + 1e-6) / (counts_nzGE['T0_M'] + 1e-6))
    LFC_sum_t = LFC_1t + LFC_2t

    hissFMt, _ = np.histogram(LFC_sum_t, bins)
    percFMt = 100 * hissFMt / hissFMt.sum() if hissFMt.sum() > 0 else np.zeros_like(hissFMt)

    hiss1t, _ = np.histogram(LFC_1t, bins)
    hiss2t, _ = np.histogram(LFC_2t, bins)
    perc1t = 100 * hiss1t / hiss1t.sum() if hiss1t.sum() > 0 else np.zeros_like(hiss1t)
    perc2t = 100 * hiss2t / hiss2t.sum() if hiss2t.sum() > 0 else np.zeros_like(hiss2t)

    T_lfc_nzGE = pd.DataFrame({
        'gRNA': T_nzGE.iloc[:, 0].values,
        'gene': T_nzGE.iloc[:, 1].values,
        'lfc1': LFC_1t.values,
        'lfc2': LFC_2t.values
    })

    gzn = [len(zGE_set), len(nzGE_set)]

    return (
        T_lfc_zGE, T_lfc_nzGE,
        bins,
        hiss1, perc1, hiss2, perc2, hissFM, percFM,
        hiss1t, perc1t, hiss2t, perc2t, hissFMt, percFMt,
        gzn
    )