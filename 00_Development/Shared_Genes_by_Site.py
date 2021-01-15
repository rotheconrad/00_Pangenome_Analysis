#!/usr/bin/env python

''' Script written to look for shared gene content in one group
not present in another group within a species. Groups are populations
at two different sites or treatments such as salinity.

The idea is to find "locally adapted gene content" or a local core
genome separate from the global core genome of species. These could be
a signal of subpopulation structure.

Input is a cd-hit cluster binary file.

Relies on the naming the scheme of genomes in the header row.
This can be adjusted on lines xxx

-------------------------------------------
Author :: Roth Conrad & Carlos Ruiz
Email :: rotheconrad@gatech.edu, cruizperez3@gatech.edu
GitHub :: https://github.com/rotheconrad, https://github.com/cruizperez/
Date Created :: July 21st, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad & Carlos Ruiz
All rights reserved
-------------------------------------------
'''

import argparse, copy
import pandas as pd


def get_subsets(df, sub, core_threshold):
    '''Collects the core genome for each subset'''
    df2 = df.filter(regex=sub)

    # set total or length of genomes (entries) in the set
    n = len(df2.columns)
    # set value to consider core
    c = n * core_threshold
    # get core gene rows
    core_df = df2[df2.sum(axis=1) >= c]
    # list of core gene names (clusters)
    core = list(core_df.index)

    if sub == '.*?': sub = 'total set'
    print(f"\nNumber of core genes in {sub} subset: {len(core)}\n")

    return core


def find_subsets(total_core, data):
    '''Finds differences in core genes between subsets'''
    
    # create container to store the unique local core genes for each subset
    unique_local_cores = {}

    # create set from the total core
    tcore = set(total_core)

    # Store the unique genes in each subset core
    subdata = {}

    # find difference(unique genes) between subset cores vs total core
    for sub, genes in data.items():
        x = list(set(genes).difference(tcore))
        subdata[sub] = x

    # find difference between subset core and other subset cores
    for sub, genes in data.items():
        x = set(genes).difference(tcore)

        # Collect genes from other subsets
        sdx = copy.deepcopy(subdata)
        sdx.pop(sub)

        subcore = []
        for k, v in sdx.items():
            subcore.extend(v)

        # these are the unique core genes to each subset
        adapted_core = list(x.difference(set(subcore)))

        # store the unique genes
        unique_local_cores[sub] = adapted_core

        
        #for gene in adapted_core: print(gene)

    return unique_local_cores


def find_other_genomes(unique_local_cores, df):
    ''' Finds other genomes outside the subset that contain the unqiue
        local core genes of the subset '''

    # create data structure to store the other genomes for each subset
    other_genomes = {}

    # get the genome names from the binary file header
    gnames = df.columns.tolist()

    # For each gene in each local subset core get other genomes with gene
    for sub, genes in unique_local_cores.items():

        # Get current subset from the f
        df2 = df.filter(regex=sub)
        # Remove the current subset from the df
        df3 = df.drop(df.filter(regex=sub).columns,axis=1)
        # define total number of genomes inside and outside the current set
        subset_total = len(df2.columns)
        other_total = len(df3.columns)

        # print header
        print(f'\n\nUnique Sub-Local Core Genes for {sub} subset:')
        print(f'Number of Genomes in {sub} subset: {subset_total}')
        print(f'Number of genomes outsite {sub} subset: {other_total}\n')
        print(
            f'Gene\tSubset_Count\tSubset_Percent\t'
            f'Other_Count\tOther_Percent'
            )
        # for each row, select genome names (columns) with the the gene (==1)
        for gene in genes:
            subset_count = len(df2.loc[gene][(df2.loc[gene] ==1)])
            subset_percent = subset_count / subset_total * 100

            other_row = df3.loc[gene][(df3.loc[gene] ==1)]
            other = other_row.index.tolist()
            other_count = len(other)
            other_percent = other_count / other_total * 100

            print(
                f'{gene}\t{subset_count}\t{subset_percent:.2f}%\t'
                f'{other_count}\t{other_percent:.2f}%'
                )



def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--binary_matrix_tsv_file',
        help='Please specify the binary matrix input tsv file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-s', '--names_to_subset',
        help='List of keywords in genome names to create subsets for!',
        metavar='',
        type=str,
        nargs='+',
        required=True
        )
    parser.add_argument(
        '-c', '--set_core_threshold',
        help='(Optional) Set the cutoff for core gene class (Default = 0.9).',
        metavar='',
        type=float,
        required=False,
        default=0.9
        )
    args=vars(parser.parse_args())

    # Read in the binary matrix tsv file as pandas dataframe
    infile = args['binary_matrix_tsv_file']
    df = pd.read_csv(infile, sep='\t', header=0, index_col=0)

    # define core threshold
    core_threshold = args['set_core_threshold']

    # store data
    data = {}

    # Collect list of core genes for entire dataset
    # This is used to compare against core genes of subsets to
    # the unique locally adapted core genes for the subset.
    total_core = get_subsets(df, '.*?', core_threshold)

    # collect core genes for each subset
    subsets = args['names_to_subset']
    for sub in subsets:
        subset = get_subsets(df, sub, core_threshold)
        data[sub] = subset

    # find differences between subset core genes
    unique_local_cores = find_subsets(total_core, data)

    # look if other genomes outside the subset have genes belonging to 
    # the local core for each subset
    _ = find_other_genomes(unique_local_cores, df)


if __name__ == "__main__":
    main()