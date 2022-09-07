"""
mutation_enrichment.py
Jorge Holguin
October 23, 2020

Set of functions (`string2range()`, `find_mutations()` and `map_mutations()`) that 
facilitate the mapping of mutations in the amino acid sequence of a protein 
(missense mutations, fusion breakpoint events) to regions of interest in a 
protein. Also includes a set of functions to calculate the enrichment of mutations
in the region (or regions) of interest compared to all other parts of a protein 
(`enrichment_analysis()` and `distribution_analysis()`). Supports the plotting of
the results with `percent_bar()` and `mean_bar()`.

Usage example:
    
    # Read the required files
    >>> df_prot = pd.read_csv(proteins_file)
    >>> df_mut = pd.read_csv(mutations_file)
    
    # Convert strings of regions to ranges or list of ranges of regions
    >>> df_prot['region_range'] = df_prot['region_string'].apply(string2range())
    
    # Find mutations in the proteins of interest
    >>> df_prot['mut_list'] = df_prot.apply(lambda x: find_mutations(x['ENST'], x['mut_list'], df_mut, identifier='Transcript_ID', 
                                            col_name = 'Mutation AA integer'), axis = 1)
    
    # Map mutations to the region (or regions) of interest (creates a new df with the results of the mapping)
    >>> df_prot_map = df_prot.apply(lambda x: map_mutations(x['mut_list'], x['region_range'], x['protein_length'], repetition = False), 
                                    axis = 1, result_type='expand')
                                    
    # Rename the columns of the new df
    >>> df_prot_map = df_prot_map.rename({0:'mut_in_region', 1:'mut_not_in_region', 2:'mut_list_in_region', 3:'mut_list_not_in_region', 
                                         4:'region_len', 5:'outside_len', 6:'prot_mut_rate', 7:'norm_region_mut_rate', 
                                         8:'norm_outside_mut_rate'}, axis = 1)
                                         
    # Concatenate df_prot and df_prot_map
    >>> df_prot = pd.concat([df_prot, df_prot_map], axis = 1)
    
    # Determine the enrichment of mutations inside the region of interest (hypergeometric test)
    >>> hypergeom_results = enrichment_analysis([df_prot], ['mut_in_region', 'mut_not_in_region', 'region_len', 'outside_len'],
                                                repetition = False)
    
    # Plot the results
    >>> percent_bar([df_prot], ['mut_in_region', 'mut_not_in_region', 'region_len', 'outside_len'], ['All proteins'], hypergeom_results,
                    ['region_hit', 'outside_hit'], verbose_labels=True, figure_size = (10, 6))
                    
    # Note that df_prot can be subset to include only specific groups of proteins. A list of dataframes can 
    # be passed onto enrichment_analysis() and percent_bar() to obtain results from subtsets of proteins

"""

import pandas as pd
from scipy.stats import hypergeom, ttest_ind, binom_test, wilcoxon, chisquare
import matplotlib.pyplot as plt

def string2range(x):
    
    """
    This function takes in a `string` representing a region of interest in a
    protein. The region of interest can be a single region or multiple regions
    of a protein. Returns a range for single regions or a list of ranges for
    multiple regions.
    
    Parameters:
    
        x (string): String containing a region or several regions of interest in a 
            protein.
            Format of x: single region -> 'start-end'
                         multiple regions -> 'start1-end1,start2-end2'
                     
    Returns:
    
        range or list of ranges: For single region proteins a range is returned. For 
            multiple region proteins a list of ranges is returned

            Format: single region -> range(start, end+1)
                    multiple region -> [range(start1, end1+1), range(start2, end2+1)]
    """
    # Handle instances with more than one range
    if ',' in x:
        list_temp = x.split(sep = ',')
        for y in range(len(list_temp)):
            list_temp[y] = list_temp[y].split(sep = '-')
        for y in range(len(list_temp)):
            for x in range(len(list_temp[y])):
                list_temp[y][x] = int(list_temp[y][x])

        # Make a range object with the bounds of the range. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        for y in range(len(list_temp)):
            for x in range(len(list_temp[y])):        
                list_temp[y] = list(range(list_temp[y][x], list_temp[y][x+1]+1))
                break

        return list(set([item for sublist in list_temp for item in sublist]))

    # Handle instances with only one range
    else:
        list_temp = x.split(sep = '-')
        for y in range(len(list_temp)):
            list_temp[y] = int(list_temp[y])

        # Make a range object with the bounds of the region. Note to the 
        # end a 1 has to be added in order to include the last position in the range
        return list(range(list_temp[0], list_temp[1]+1))
    

def find_mutations(ids, mut_list, df_mut, identifier, col_name_mut, col_name_rec):
    
    """
    Given an identifier `ids` find the mutations associated with 
    that identifier from `df_mut` and return the list of mutations found
    plus the already present mutations (`mut_list`). If no new mutations are found
    then return the already present mutations (`mut_list`).
    
    Parameters:
    
        ids (string): Identifier for the protein of interest. The identifier used
            will depend on the identifiers available in `df_mut`. Can typically be UniProt,
            ENSP or ENST. The use of Ensembl identifiers over Uniprot identifiers is 
            recommended to avoid retrieving mutations from multiple isoforms of a protein.

        mut_list (list): List of existing mutations for the protein of interest

        df_mut (pandas dataframe): Data frame containing mutations. Must contain
            a column with the same name as `identifier`

        identifier (string): Name of the column in `df_mut` that contains identifiers.
            Must be of the same format as the `ids` for the protein of interest.
            
        col_name (string): Name of the column in `df_mut` that contains the mutations
        
    Returns:
    
        list: If new mutations are found, returns a list of new mutations + old mutations 
            If no new mutations are found, returns the list of old mutations
    """
    
    df_temp = df_mut.loc[df_mut[identifier] == ids]
    df_temp['Mut_aa,Mut_recurrency'] = list(zip(df_temp[col_name_mut], df_temp[col_name_rec]))              

    temp_list = list(df_temp['Mut_aa,Mut_recurrency'])

    if len(temp_list) > 0:
        
        return temp_list + mut_list
        
    else:
        
        return mut_list

    
def map_mutations(mut_list, region_range, prot_len, repetition, region_type = 'range'):
    
    """
    Maps the mutations given in `mut_list` to the region of a protein given in 
    `region_range` and to outside the region of interest. The variable `repetition` is
    boolean and indicates whether all the mutations are to be considered (`repetition` == True)
    or if only the unique mutations are to be considered (`repetition` == False).
    
    Parameters:
    
        mut_list (list): List of mutations for a given protein

        region_range (range or list): Range of the region of interest range(start, end) or
            list of ranges with regions of interest [range(start1, end1), range(start2, end2)]

        prot_len (int): Length of the protein in amino acids

        repetition (bool): If True all mutations are considered when mapping. If False
            only unique mutations are considered when mapping.
        
    Returns:
    
        list: Contains the calculated values explained below

            hit_count (int): Number of mutations inside the region of interest

            miss_count (int): Number of mutations outside the region of interest

            mut_in_region (list): List of mutations inside the region of interest

            mut_not_in_region (list): List of mutations outside the region of interest

            len(region_res) (int): Number of residues in the region of interest

            len(outside_res) (int): Number of residues outside the region of interest

            prot_mut_rate (float): Mutation rate of the protein (number of mutations/protein length)
            
            region_mut_rate (float): Mutation rate of the region of interest (number of mutations/length of region)
            
            outside_mut_rate (float): Mutation rate outside the region of interest (number of mutations/length outside region of interest)

            norm_region_mut_rate (float): Normalized mutation rate in the region of interest
                (hit_count/len(region_res))/prot_mut_rate

            norm_outside_mut_rate (float): Normalized mutation rate outside the region of interest
                (miss_count/len(outside_res))/prot_mut_rate
    """
    
    if region_type == 'range':
    
        # Here the proteins with only one IAS-range are analyzed
        if type(region_range) == range:

            # Create a set with the residues inside and outside the IAS-range
            region_res = set(region_range)
            protein_res = set(range(1, prot_len + 1))
            outside_res = protein_res - region_res

        # Here proteins with more than one IAS are analyzed    
        elif type(region_range) == list:

            # Set a variable with the list of IAS ranges
            list_region_range = region_range

            # All the residues in all the IAS ranges for a protein are added into
            # one list
            region_res = []

            for item in range(len(list_region_range)):
                region_res = region_res + list(list_region_range[item])

            # Create a set with the residues outside the IAS range
            region_res = set(region_res)
            protein_res = set(range(1, prot_len + 1))
            outside_res = protein_res - region_res
            
    elif region_type == 'interface':
        
        if type(region_range) == list:

            # Create a set with the residues outside the IAS range
            region_res = set(region_range)
            protein_res = set(range(1, int(prot_len) + 1))
            outside_res = protein_res - region_res
        
    # Get the list of mutations for a protein
    if repetition == True:
        mutation_list = mut_list
    elif repetition == False:
        # Here we set the recurrency of the mutation to 1 (no repetition)
        mutation_list = [(i[0], 1) for i in mut_list]
    
    # Count the mutations that map inside and outside the IAS range
    hit_count = 0
    miss_count = 0
    
    # Make two lists where the mutations inside and outside the IAS are stored
    mut_in_region = []
    mut_not_in_region = []
    
    # Iterate through the mutation list and determine if it occurs inside
    # or outside the IAS
    for mut in mutation_list:
        if mut[0] in region_res:
            hit_count = hit_count + mut[1]
            mut_in_region.append(mut)
            
        elif mut[0] in outside_res:
            miss_count = miss_count + mut[1]
            mut_not_in_region.append(mut)
    
    p_value = binom_test(hit_count, hit_count + miss_count, len(region_res)/len(protein_res), alternative = 'greater')
    
    prot_mut_rate = (hit_count + miss_count)/prot_len
    
    if prot_mut_rate > 0:
        
        if len(region_res) > 0:
            region_mut_rate = hit_count/len(region_res)
            norm_region_mut_rate = (hit_count/len(region_res))/prot_mut_rate
            
        else: 
            region_mut_rate = 0
            norm_region_mut_rate = 0
        
        if len(outside_res) > 0:
            outside_mut_rate = miss_count/len(outside_res)
            norm_outside_mut_rate = (miss_count/len(outside_res))/prot_mut_rate
            
        else:
            outside_mut_rate = 0
            norm_outside_mut_rate = 0
        
    else: 
        
        region_mut_rate = 0
        outside_mut_rate = 0
        
        norm_region_mut_rate = 0 
        norm_outside_mut_rate = 0 
    
    return [hit_count, miss_count, mut_in_region, mut_not_in_region, len(region_res), len(outside_res), prot_mut_rate, region_mut_rate,
            outside_mut_rate, norm_region_mut_rate, norm_outside_mut_rate, p_value]


def enrichment_analysis(dfs, cols, repetition):
    
    """
    Takes a list of dataframes `dfs` and the boolean parameter `repetition` that indicates if the 
    enrichment analysis is to be done allowing for repetition (binomial) or not allowing for repetition
    (hypergeometric). Also takes in a list of columns `cols` from the `dfs` that are used to compute the
    enrichment analysis.
     
    Parameters:
    
        dfs (list): List of dataframes to be used in the statistical analysis

        cols (list): List of strings with the name of the columns in the `dfs` to be used for statistical analysis
            Format of cols: [mutations in region, mutations outside region, length of region, length outside region]

        repetition (bool): Determines if the enrichment analyis is done allowing for repetition (binomial)
            or not (hypergeometric)
        
    Returns:
    
        tuple: List of p-values (one per df in dfs) and a list of annotations that can be used for plotting

            Format of annotations: 

                NS -> not significat at a 0.05 threshold
                '*' -> 0.05 > p-value > 0.01
                '**' -> 0.01 > p-value > 0.001
                '***' -> 0.001 > p-value > 0.0 
    """
    
    p_values = []
    annot_text = []
    
    # In the hypergeometric test, k = number of successes, M = size of the population, 
    # n = number of successes in the population, N = sample size (draws)
    if repetition == False:
        
        for df in dfs:
            
            proportion_mut_region = sum(df[cols[0]])/sum(df[cols[2]])
            proportion_mut_outside = sum(df[cols[1]])/sum(df[cols[3]])
            
            if proportion_mut_region > proportion_mut_outside:
                
                p_values.append(hypergeom.sf(sum(df[cols[0]])-1, sum(df[cols[3]])+sum(df[cols[2]]), sum(df[cols[2]]), sum(df[cols[0]]) + sum(df[cols[1]])))
                
            elif proportion_mut_outside > proportion_mut_region:
                
                p_values.append(hypergeom.cdf(sum(df[cols[0]]), sum(df[cols[3]])+sum(df[cols[2]]), sum(df[cols[2]]), sum(df[cols[0]]) + sum(df[cols[1]])))
     
    
    # Binomial test (successes, trials, p(success))
    elif repetition == True:
        
        for df in dfs:
            
            p_values.append(binom_test(sum(df[cols[0]]), sum(df[cols[0]]) + sum(df[cols[1]]), sum(df[cols[2]])/(sum(df[cols[3]])+sum(df[cols[2]]))))
         
    for value in p_values:
        if value > 0.05:
            annot_text.append('NS')
        elif value <= 0.05 and value > 0.01:
            annot_text.append('*')
        elif value <= 0.01 and value > 0.001:
            annot_text.append('**')
        elif value <= 0.001 and value > 0:
            annot_text.append('***')
        elif value == 0.0:
            annot_text.append('***')
    
    return (p_values, annot_text)
        

def percent_bar(dfs, cols, x_labels, stat_result, legend, verbose_labels, figure_size, save_path = '', barwidth = 0.1, 
                colors = ['#ff0000', '#808080'], annotate = 'stars'):
    
    """
    Plots grouped percentage bars indicating the percent number of mutated residues 
    inside a region of a protein and the percent number of mutated residues outside 
    that region. This is done to determine whether there is an enrichment of mutations
    inside the region of interest. Supports the plotting of data from multiple pandas
    data frames corresponding to distinct groups of proteins.
    
    Parameters:
    
        dfs (list): List of dataframes to use for plotting

        cols (list): List of strings indicating the column names used in plotting. Must follow the format below
            Format of cols: [Mutations in region, Mutations outside region, Region Length, Outside Length]

        x_labels (list): List of strings indicating the names of the groups of proteins. Must be the same
            length as `dfs`

        stat_result (tuple): Tuple containing a list of p-values from the statistical analyis and a list of the 
            graphical annotations to use when plottig. ie ([0.06, 0.05], ['NS', '*']).

        legend (list): List of strings indicating the names of the legends used when plotting. Must be of
            length 2. ie ['legend1', 'legend2']

        verbose_labels (bool): Indicates whether to include verbose labels in the x-axis. If True the labels
            in the x-axis will contain sample size of genes, sample size of mutations and p-value. If False
            only the names in `x_labels` will be used.

        figure_size (tuple): Size of the figure (x-size, y-size).

        save_path (string): Path where to save the figure. The figure is not saved by default

        barwidth (float): Size of the bars, 0.1 by default

        colors (list): Color identifiers to use when plotting. Must be of length 2. Red and grey are the default

        annotate (string): Either 'stars' (default) or 'p_value'. Controls the label on top of the error bars.
            Uses the information inside `stat_results`.
    
    Returns:
    
        matplotlib.pyplot.plot: Returns the plot.
    
    """
    
    bars1 = []
    bars2 = []
    sample_size_genes = []
    sample_size_mutations = []
    value_pairs = []
    tick_labels = []
    
    r1 = []
    
    i = 0
    
    for df in dfs:
        
        percent_in = (sum(df[cols[0]])/sum(df[cols[2]]))*100
        percent_out = (sum(df[cols[1]])/sum(df[cols[3]]))*100
        number_genes = len(df)
        number_mutations = sum(df[cols[0]]) + sum(df[cols[1]])
        
        bars1.append(percent_in)
        bars2.append(percent_out)
        sample_size_genes.append(number_genes)
        sample_size_mutations.append(number_mutations)
        value_pairs.append((percent_in, percent_out))
        
        r1.append(barwidth*(2.5*i))
        
        if verbose_labels == True:
            
            tick_labels.append(x_labels[i] + '\nproteins = %d\nmutations = %d\np-value = %s' % (number_genes, number_mutations, 
                                                                                                      "{:.2e}".format(stat_result[0][i])))
            
        i = i + 1
    
    r2 = [x + barwidth for x in r1]
    
    plt.figure(figsize=figure_size)
    
    plt.bar(r1, bars1, color=colors[0], width=barwidth, edgecolor='white', label=legend[0])
    plt.bar(r2, bars2, color=colors[1], width=barwidth, edgecolor='white', label=legend[1])
    plt.ylabel('Percentage')
    
    if verbose_labels == True:
        plt.xticks([r + barwidth/2 for r in r1], tick_labels)
    else:
        plt.xticks([r + barwidth/2 for r in r1], x_labels)
    
    for i in range(len(r1)):
        x1, x2 = r1[i], r2[i]
        y, h, col = max(value_pairs[i]) + max(value_pairs[0])*0.04, max(value_pairs[0])*0.04, 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.0, c=col)
        
        if annotate == 'stars':
            plt.text((x1+x2)*0.5, y+h, stat_result[1][i], ha='center', va='bottom', color=col)
        elif annotate == 'p_value':
            plt.text((x1+x2)*0.5, y+h, "{:.2e}".format(stat_result[0][i]), ha='center', va='bottom', color=col)
            
#     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.legend(loc='upper left')
    plt.tight_layout()
    
    if len(save_path) > 0:
        plt.savefig(save_path, bbox_inches='tight', dpi = 300)
        
     
def enrichment_bar(df, index, column, figure_size, save_path = ''):

    """
    Plots grouped percentage bars indicating the percent number of mutated residues 
    inside a region of a protein and the percent number of mutated residues outside 
    that region. This is done to determine whether there is an enrichment of mutations
    inside the region of interest. Supports the plotting of data from multiple pandas
    data frames corresponding to distinct groups of proteins.
    
    Parameters:
    
        df (pandas dataframe): Dataframes to use for plotting
        
        index (string): Name of the column to use for the x-axis
        
        column (string): Name of the column to use for the y-axis

        figure_size (tuple): Size of the figure (x-size, y-size).

        save_path (string): Path where to save the figure. The figure is not saved by default
    
    Returns:
    
        matplotlib.pyplot.plot: Returns the plot.
    
    """    

    # Make a df where the index are the gene names and the colums are the cancer types
    pivot_df = df.pivot(index = index, columns = column, values = column)
    
    # Make a column with the sum of all colums for a given gene (helps with sorting)
    pivot_df['Sum'] = pivot_df.sum(axis = 1, skipna = True)
    
    # Sort the df from highest to lowest number in the Sum column
    pivot_df = pivot_df.sort_values(by = 'Sum', ascending=False)
    
    # Drop the Sum column
    pivot_df = pivot_df.drop(['Sum'], axis=1)
    
    # Plot the results
    pivot_df[:40].dropna(axis='columns', how = 'all').plot.bar(stacked=True, color = 'k', edgecolor = 'k', 
                                                               width = 0.8, figsize=figure_size, legend = None)
    plt.ylabel(column)
    plt.xlabel('Gene')
    
    if len(save_path) > 0:
        plt.savefig(save_path, bbox_inches='tight', dpi = 300)