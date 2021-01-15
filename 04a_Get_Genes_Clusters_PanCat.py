#!/usr/bin/env python

'''Generate .tsv file Genes Clusters Pangenome Categories

Takes the Binary Matrix .tsv file from 03 and the CD-HIT cluster file 
input, and outputs a .tsv file of:
Gene_Name Cluster_Name Pangenome_Category
Where Pangenome_Category is one of: (Core, Common, Rare, Specific)

This tool takes the following input parameters:

    * bnry - binary gene presence absence matrix .tsv file
    * clstr - CD-HIT cluster file
    * out - output file name

This script returns the following files:

    * outputfilename.tsv

This script requires the following packages:

    * argparse

This file can also be imported as a module and contains the follwing 
functions:

    * get_category - reads the bnry file and assigns pangenome category
    * build_the_list - Orchestrates Parsing, computing, and output
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, September 5th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


def get_category(bnry):
    ''' Read bnry file and calculate the category of each cluster '''

    pan_category = {} # by cluster. Core, Common, Rare, Specific 

    with open(bnry, 'r') as f:
        header = f.readline().rstrip().split()
        genomes = header # genome names are the header
        n = len(genomes) # This is the number of genomes ie(100)
        # Define gene category cutoffs
        core, rare, specific = 0.90, 0.20, 1/n
        # each cluster is a line in the bnry file
        # compute gene category for each cluster
        for l in f:
            X = l.rstrip().split('\t')
            cluster = X[0] # Cluster # is the first column of bnry
            gs = [int(i) for i in X[1:]] # genes in cluster count
            g = sum(gs) # This is how many genomes cluster is in
            v = g / n # percent genomes with gene
            if v == specific: pan_category[cluster] = ['Specific', v]
            elif v < rare: pan_category[cluster] = ['Rare', v]
            elif v < core: pan_category[cluster] = ['Common', v]
            elif v >= core: pan_category[cluster] = ['Core', v]
            else: print(l, v, 'error in assigning gene category')

    return pan_category

def build_the_list(bnry, clstr, out):
    ''' Reads in the files, builds the list, and writes to out '''

    pan_category = get_category(bnry)

    with open(clstr, 'r') as f, open(out, 'w') as o:
        o.write('Gene_Name\tCluster_Name\tPangenome_Category\tn/N\n')
        # n/N = number of genomes with gene in cluster / total genomes

        for l in f:
            if l.startswith('>'): cluster = '_'.join(l.rstrip()[1:].split(' '))
            else:
                gene = l.split('>')[1].split('...')[0]
                pcat = pan_category[cluster]
                line_out = f'{gene}\t{cluster}\t{pcat[0]}\t{pcat[1]}\n'
                o.write(line_out)


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-bnry', '--binary_file',
        help='Please specify the binary.tsv input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-clstr', '--cluster_file',
        help='Please specify the CD-HIT Cluster File!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-out', '--output_file',
        help='Please specify the name to us for the output .tsv file',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Generating Genes, Clusters, and PanCat list...')

    build_the_list(
        args['binary_file'],
        args['cluster_file'],
        args['output_file']
        )

if __name__ == "__main__":
    main()
