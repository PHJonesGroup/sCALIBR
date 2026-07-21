import numpy as np
import pandas as pd

def general_volcano(alf, sfdr_corr, thr_scoreh, thr_scored, score, fdr, cond, genes, plot=True, ax=None):
    thr_fdr = alf
    fdr_corr = fdr + sfdr_corr  # Avoid log10(0)
    LPV = -np.log10(fdr_corr)

    T_gene = pd.DataFrame({'gene': genes, 'score': score, 'fdr_corr': fdr_corr})

    indh = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s >= thr_scoreh)]
    indd = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s <= thr_scored)]

    T_hs = T_gene.loc[indh].copy() if indh else pd.DataFrame()
    T_ds = T_gene.loc[indd].copy() if indd else pd.DataFrame()

    if plot and ax is not None:
        ma = max(5, np.max(LPV) + 1 if LPV.size > 0 else 5)
        ax.plot(score, LPV, 'pk', markersize=4, label='All genes')
        ax.grid(True)
        ax.set_ylabel('-Log10 FDR', fontsize=14)
        ax.set_xlabel('score', fontsize=14)
        ax.set_ylim([0, ma])

        if not T_hs.empty:
            ax.plot(T_hs['score'], -np.log10(T_hs['fdr_corr']), 'pr', linewidth=2, label='Hits')
            for i, gene_name in enumerate(T_hs['gene'][:20]):
                ax.text(T_hs['score'].iloc[i], -np.log10(T_hs['fdr_corr']).iloc[i], gene_name, va='bottom', ha='right')

        if not T_ds.empty:
            ax.plot(T_ds['score'], -np.log10(T_ds['fdr_corr']), 'pc', linewidth=2, label='Depleted')
            for i, gene_name in enumerate(T_ds['gene'][:20]):
                ax.text(T_ds['score'].iloc[i], -np.log10(T_ds['fdr_corr']).iloc[i], gene_name, va='bottom', ha='right')

        ax.legend()

    return LPV, indh, indd, T_ds, T_hs, T_gene
