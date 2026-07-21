def count_vertical_names(gene_names):
    """
    Given a list of gene names (strings), this function groups consecutive identical names,
    returning:
    - ggenes: list of unique gene names for each consecutive group
    - gg: list of counts of repeats per unique gene name group
    - ind_gn: list of start indices (1-based) for each group in the original list
    
    Example:
    gene_names = ['Aars2', 'Aars2', 'Aars2', 'Aars2', 'Aasdhppt', 'Aasdhppt', 'Aasdhppt', 'Aasdhppt', ...]
    ggenes = ['Aars2', 'Aasdhppt', ...]
    gg = [4, 4, ...]  # counts
    ind_gn = [1, 5, ...]  
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
