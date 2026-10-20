import pandas as pd

def summarise(T_vert, colnames, norm_method):
    """Return per-gene median, mean, std for one quantity."""
    long = T_vert.melt(id_vars=['gene'], value_vars=colnames,
                        var_name='rep', value_name='v')
    long['v'] = pd.to_numeric(long['v'], errors='coerce')

    if norm_method == 'pooled':
        flat = long                                   # all guide x rep values
    elif norm_method == 'sample':
        # collapse guides within each replicate first
        flat = (long.groupby(['gene', 'rep'], as_index=False)['v']
                    .median())
    else:
        raise ValueError(f"norm_method must be 'pooled' or 'sample', "
                            f"got '{norm_method}'")

    g = flat.groupby('gene')['v']
    return g.median(), g.mean(), g.std(ddof=1)