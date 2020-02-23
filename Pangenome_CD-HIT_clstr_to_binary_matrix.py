#!/usr/bin/env python

''' Parse CD-HIT clstr into binary gene presence/absence matrix.

This script parses a CD-HIT clstr file into a binary matrix of gene 
cluster presence/absence. Genomes are written as the columns and gene
clusters are written as the rows.

This script is written for clustered genes from Prodigal gene predictions
where fasta sequence names are of the form >UniqueGeneID_Contig#_Gene#.
The script expects underscore delimited sequence names where when split,
positionally UniqueGeneID is in position -3, Contig# -2 and Gene# -1.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 24th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import pandas as pd

def generate_binary_matrix(GenomeIDs, ClusterIDs):

    print('Generating binary matrix ...')

    cid = list(ClusterIDs.keys())

    matrix = defaultdict(list)

    for genome, clusters in GenomeIDs.items():
        for c in cid:
            if c in clusters: matrix[genome].append(1)
            else: matrix[genome].append(0)

    binary_matrix = pd.DataFrame(matrix, index=cid)

    return binary_matrix


def parse_cdhit_clstr(clstr):

    print('Parsing CD-HIT Clstr File ...')
    # Initialize variables.
    ClusterIDs = defaultdict(lambda: defaultdict(int))
    GenomeIDs = defaultdict(lambda: defaultdict(int))

    clusterID = None

    # Read through clstr file and populate initialized variables.
    with open(clstr, 'r') as c:
        for l in c:
            if l.startswith('>'): # This denotes the beginning of a cluster.
                if clusterID: # When cluster is set, populate dictionary.
                    clusterID = None # Reset for new cluster
                clusterID = l[1:8] + '_' + l[9:].rstrip() # name new cluster.

            elif l[0].isdigit(): # This denotes an entry within a cluster
                X = l.rstrip().split(' ') # Split line by space

                # Select Unique Genome ID
                genomeID = X[1].split('_')[-3]

                GenomeIDs[genomeID][clusterID] += 1
                ClusterIDs[clusterID][genomeID] += 1

    return GenomeIDs, ClusterIDs

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-c', '--cdhit_clstr_file',
        help='Please specify a CD-HIT clstr file for input!',
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='What would you like to call the output tsv file?',
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')

    GenomeIDs, ClusterIDs = parse_cdhit_clstr(
                                            args['cdhit_clstr_file'],
                                            )

    binary_matrix = generate_binary_matrix(GenomeIDs, ClusterIDs)

    binary_matrix.to_csv(args['output_file'], sep='\t')

    print(f'Success!! Binary matrix written to file.')


if __name__ == "__main__":
    main()