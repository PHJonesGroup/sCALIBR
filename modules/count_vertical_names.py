def count_vertical_names(gene_names):
    """
    Given a list of gene names (strings), this function groups rows with identical names.
    
    Parameters
    ----------
    gene_names : list of str
        Gene names found in column gene

    Returns
    -------
    ggenes : list of str
        List of unique gene names
    gg : list of int
        Number of times each gene identifier appears in dataframe
    ind_gn : list of int
        List of start indices (1-based) for each group in the original list
    """
    ggenes = []
    gg = []
    ind_gn = []

    j = 0  
    k = 1  
    
    ggenes.append(gene_names[0])
    ind_gn.append(1) 
    
    for i in range(len(gene_names) - 1):
        s1 = gene_names[i]
        s2 = gene_names[i + 1]
        
        if s1 == s2:
            k += 1
        else:
            gg.append(k)
            ggenes[j] = gene_names[i]  
            j += 1
            k = 1
            ggenes.append(gene_names[i + 1]) 
            ind_gn.append(i + 2)  
    
    gg.append(k)
    ggenes[j] = gene_names[-1]

    return ggenes, gg, ind_gn
