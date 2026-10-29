import numpy as np
import pandas as pd
from scipy.special import gammaincc

from .bh_adjust import bh_adjust


def fisher_combine(p, scale_factor=1/3):
    """
    Combine p-values with Fisher's method, scaled for possible dependency.

    chi = -2 * sum(log p) * scale_factor, df = 2 * n * scale_factor, and the
    combined p is the chi-square upper tail, gammaincc(df / 2, chi / 2).
    p-values of exactly 0 are set to machine epsilon; NaN values are ignored.

    Parameters
    ----------
    p : array-like
        p-values to combine (e.g. one gene's gRNAs).
    scale_factor : float, optional
        Down-weights the evidence when the p-values may be dependent; 1 is
        plain Fisher, 1/3 a rough allowance for dependency.

    Returns
    -------
    chi_stat : float
        Scaled chi-square statistic (NaN if no finite p-values).
    combined_p : float
        Combined p-value (NaN if no finite p-values).
    """
    p = np.asarray(p, dtype=float).ravel()
    p = p[np.isfinite(p)]
    if p.size == 0:
        return np.nan, np.nan
    p = np.where(p == 0, np.finfo(float).eps, p)
    chi_stat = -2 * np.sum(np.log(p)) * scale_factor
    df = 2 * p.size * scale_factor
    return chi_stat, gammaincc(df / 2, chi_stat / 2)


def _rep_suffixes(T_vert):
    """Replicate suffixes shared by the normalised_lfc and Q columns ('' or '_<rep>')."""
    sfx = [c[len('Q'):] for c in T_vert.columns if c == 'Q' or c.startswith('Q_')]
    if not sfx:
        raise KeyError(f"No Q columns. Columns: {T_vert.columns.tolist()}")
    return sfx


def _per_gene(T_vert, nl_cols, p_cols, scale_factor):
    """Median normalised LFC and Fisher-combined p per gene over the given columns."""
    rows = []
    for gene, g in T_vert.groupby('gene', sort=True):
        nl = g[nl_cols].to_numpy(dtype=float).ravel()
        p = g[p_cols].to_numpy(dtype=float).ravel()
        _, fp = fisher_combine(p, scale_factor)
        rows.append({'gene': gene,
                     'normalised_lfc': np.nanmedian(nl) if np.isfinite(nl).any() else np.nan,
                     'fisher_p': fp,
                     'n_p': int(np.isfinite(p).sum())})
    return pd.DataFrame(rows)


def gene_fisher(T_vert, norm_method="pooled", scale_factor=1/3):
    """
    Gene-level normalised LFC and Fisher-combined p-values, BH-corrected across genes.

    'pooled' : per gene, the median normalised LFC over all guide x replicate
               values, and Fisher's method over all their gRNA p-values (Q);
               BH over the Fisher p-values of all target genes.
    'sample' : per replicate, the median normalised LFC and Fisher's method
               over the gene's guides in that replicate, BH within that
               replicate (one result per sample). Then a final result pooling
               the samples: Fisher's method once over all the gene's raw gRNA
               p-values from every sample (guides x samples), BH across genes,
               and the median of the per-replicate median normalised LFCs.
               Fisher's method is never applied to Fisher p-values.

    Parameters
    ----------
    T_vert : pandas.DataFrame
        Per-gRNA table [gRNA, gene, gene_type, lfc(_<rep>)...,
        normalised_lfc(_<rep>)..., Q(_<rep>)...], targets only.
    norm_method : {'pooled', 'sample'}
        Aggregation mode, as above.
    scale_factor : float
        Passed to fisher_combine.

    Returns
    -------
    pandas.DataFrame
        'pooled': [gene, normalised_lfc, fisher_p, q, n_p].
        'sample': [gene, normalised_lfc, fisher_p, q, n_p, then per replicate
        normalised_lfc_<rep>, fisher_p_<rep>, q_<rep>]; fisher_p / q are the
        final result over all samples' raw gRNA p-values.
    """
    sfx = _rep_suffixes(T_vert)

    if norm_method == "pooled":
        out = _per_gene(T_vert, [f'normalised_lfc{s}' for s in sfx], [f'Q{s}' for s in sfx],
                        scale_factor)
        out['q'] = bh_adjust(out['fisher_p'])
        return out[['gene', 'normalised_lfc', 'fisher_p', 'q', 'n_p']]

    if norm_method != "sample":
        raise ValueError(f"norm_method must be 'pooled' or 'sample', got '{norm_method}'")

    per_rep = []
    for s in sfx:
        r = _per_gene(T_vert, [f'normalised_lfc{s}'], [f'Q{s}'], scale_factor)
        r[f'q{s}'] = bh_adjust(r['fisher_p'])
        per_rep.append(r.rename(columns={'normalised_lfc': f'normalised_lfc{s}',
                                         'fisher_p': f'fisher_p{s}'})
                        .drop(columns='n_p').set_index('gene'))
    wide = pd.concat(per_rep, axis=1)

    # final result: Fisher once over all raw gRNA p-values of the gene, every sample
    out = _per_gene(T_vert, [f'normalised_lfc{s}' for s in sfx], [f'Q{s}' for s in sfx],
                    scale_factor)
    out['normalised_lfc'] = wide[[f'normalised_lfc{s}' for s in sfx]].median(axis=1).reindex(out['gene']).to_numpy()
    out['q'] = bh_adjust(out['fisher_p'])
    out = out[['gene', 'normalised_lfc', 'fisher_p', 'q', 'n_p']]
    rep_cols = [c for s in sfx for c in (f'normalised_lfc{s}', f'fisher_p{s}', f'q{s}')]
    return out.join(wide[rep_cols], on='gene')
