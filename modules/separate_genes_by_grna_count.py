import numpy as np

def separate_genes_by_grna_count(T_vert, ggenes, gene_names, ind_gn, gg, rep, n_cond=2):
    wt = T_vert.iloc[:, 2:(2 + n_cond*rep)].to_numpy()

    groups = {}   # count -> dict with genes, gene_names, indices, data
    for i in range(len(gg)):
        c = gg[i]
        start = ind_gn[i] - 1
        g = groups.setdefault(c, {"genes": [], "gene_names": [], "ind": [], "wt": []})
        g["genes"].append(ggenes[i])
        g["gene_names"].extend(gene_names[start:start + c])
        g["ind"].extend(range(start, start + c))
        g["wt"].append(wt[start:start + c, :])

    for c, g in groups.items():
        g["wt"] = np.vstack(g["wt"]) if g["wt"] else np.empty((0, wt.shape[1]))

    nums = {c: len(g["genes"]) for c, g in groups.items()}
    return groups, nums