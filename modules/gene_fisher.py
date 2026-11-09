import numpy as np
import pandas as pd
from scipy.special import gammaincc

from .bh_adjust import bh_adjust


def fisher_combine(p, scale_factor=1.0):
    """
    Combine p-values with Fisher's method, optionally scaled for dependency.

    chi = -2 * sum(log p) * scale_factor, df = 2 * n * scale_factor, and the
    combined p is the chi-square upper tail, gammaincc(df / 2, chi / 2).
    p-values of exactly 0 are set to machine epsilon; NaN values are ignored.

    Parameters
    ----------
    p : array-like
        p-values to combine (e.g. one gene's gRNAs), each running 0 to 1
        under the null.
    scale_factor : float, optional
        Down-weights the evidence when the p-values are dependent; 1 is
        plain Fisher (the default: control pseudo-genes are calibrated with it).

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
    """Replicate suffixes shared by the normalised_lfc and p columns ('' or '_<rep>')."""
    sfx = [c[len('dep_p'):] for c in T_vert.columns if c.startswith('dep_p')]
    if not sfx:
        raise KeyError(f"No dep_p columns. Columns: {T_vert.columns.tolist()}")
    return sfx


def _per_gene(T_vert, nl_cols, sfx, scale_factor):
    """
    Median normalised LFC and directional Fisher p per gene over the replicates in sfx.

    Fisher's method runs twice: on every gRNA's depletion p (is the gene
    depleted?) and on every gRNA's enrichment p (is it enriched?). A gRNA moving
    the other way has a p near 1 and counts against the gene. The gene p is
    twice the smaller of the two, capped at 1, for testing both directions.
    """
    dep_cols = [f'dep_p{s}' for s in sfx]
    enr_cols = [f'enr_p{s}' for s in sfx]
    rows = []
    for gene, g in T_vert.groupby('gene', sort=True):
        nl = g[nl_cols].to_numpy(dtype=float).ravel()
        dep = g[dep_cols].to_numpy(dtype=float).ravel()
        _, fp_dep = fisher_combine(dep, scale_factor)
        _, fp_enr = fisher_combine(g[enr_cols].to_numpy(dtype=float).ravel(), scale_factor)
        rows.append({'gene': gene,
                     'normalised_lfc': np.nanmedian(nl) if np.isfinite(nl).any() else np.nan,
                     'fisher_p': min(1.0, 2 * min(fp_dep, fp_enr)),
                     'direction': 'depleted' if fp_dep <= fp_enr else 'enriched',
                     'fisher_p_dep': fp_dep,
                     'fisher_p_enr': fp_enr,
                     'n_p': int(np.isfinite(dep).sum())})
    return pd.DataFrame(rows)


def gene_fisher(T_vert, norm_method="pooled", scale_factor=1.0):
    """
    Gene-level normalised LFC and directional Fisher p-values, BH-corrected across genes.

    For each gene, Fisher's method is run separately on the depletion and the
    enrichment p-values of its gRNAs; the gene p is twice the smaller of the
    two (capped at 1) and 'direction' records which one it was.

    'pooled' : per gene, the median normalised LFC over all guide x replicate
               values, and the directional Fisher test over all their gRNA
               p-values; BH over the gene p-values of all target genes.
    'sample' : per replicate, the median normalised LFC and the directional
               Fisher test over the gene's guides in that replicate, BH within
               that replicate (one result per sample). Then a final result
               pooling the samples: the directional Fisher test once over all
               the gene's raw gRNA p-values from every sample (guides x
               samples), BH across genes, and the median of the per-replicate
               median normalised LFCs. Fisher's method is never applied to
               Fisher p-values.

    Parameters
    ----------
    T_vert : pandas.DataFrame
        Per-gRNA table [gRNA, gene, gene_type, lfc(_<rep>)...,
        normalised_lfc(_<rep>)..., p(_<rep>), dep_p(_<rep>), enr_p(_<rep>)...],
        targets only; dep_p / enr_p = raw empirical one-sided gRNA p-values.
    norm_method : {'pooled', 'sample'}
        Aggregation mode, as above.
    scale_factor : float
        Passed to fisher_combine.

    Returns
    -------
    pandas.DataFrame
        'pooled': [gene, normalised_lfc, fisher_p, direction, q, n_p,
        fisher_p_dep, fisher_p_enr].
        'sample': the same, then per replicate normalised_lfc_<rep>,
        fisher_p_<rep>, direction_<rep>, q_<rep>; fisher_p / q are the final
        result over all samples' raw gRNA p-values.
    """
    sfx = _rep_suffixes(T_vert)
    cols = ['gene', 'normalised_lfc', 'fisher_p', 'direction', 'q', 'n_p',
            'fisher_p_dep', 'fisher_p_enr']

    if norm_method == "pooled":
        out = _per_gene(T_vert, [f'normalised_lfc{s}' for s in sfx], sfx, scale_factor)
        out['q'] = bh_adjust(out['fisher_p'])
        return out[cols]

    if norm_method != "sample":
        raise ValueError(f"norm_method must be 'pooled' or 'sample', got '{norm_method}'")

    per_rep = []
    for s in sfx:
        r = _per_gene(T_vert, [f'normalised_lfc{s}'], [s], scale_factor)
        r[f'q{s}'] = bh_adjust(r['fisher_p'])
        per_rep.append(r.rename(columns={'normalised_lfc': f'normalised_lfc{s}',
                                         'fisher_p': f'fisher_p{s}',
                                         'direction': f'direction{s}'})
                        [['gene', f'normalised_lfc{s}', f'fisher_p{s}', f'direction{s}', f'q{s}']]
                        .set_index('gene'))
    wide = pd.concat(per_rep, axis=1)

    # final result: Fisher once over all raw gRNA p-values of the gene, every sample
    out = _per_gene(T_vert, [f'normalised_lfc{s}' for s in sfx], sfx, scale_factor)
    out['normalised_lfc'] = wide[[f'normalised_lfc{s}' for s in sfx]].median(axis=1).reindex(out['gene']).to_numpy()
    out['q'] = bh_adjust(out['fisher_p'])
    out = out[cols]
    rep_cols = [c for s in sfx
                for c in (f'normalised_lfc{s}', f'fisher_p{s}', f'direction{s}', f'q{s}')]
    return out.join(wide[rep_cols], on='gene')
