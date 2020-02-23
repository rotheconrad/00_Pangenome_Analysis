#!/usr/bin/env python

'''Get Fasta Sequence for CD-HIT Cluster

Takes the *.fnn.clstr file from CD-HIT(-EST) and the concatenated
genes file used as input to CD-HIT and returns a fasta file of
sequences for the user specified gene clusters.

This tool takes the following input parameters:

    * clstr - CD-HIT cluster file (str)
    * fasta - CD-HIT input fasta file (str)
    * cnumb - Cluster number to get fasta sequences for. (int)
    * out - name for the output fasta file(str)

This script returns the following files:

    * fasta file of cnumb sequneces with name {out}

This script requires the following packages:

    * argparse

This file can also be imported as a module and contains the follwing 
functions:

    * get_seqs_in_clstr - reads CD-HIT clstr file and retrieves seq names
    * get_seq_fasta - reads fasta file and gets sequences
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: November 4th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def get_seq_fasta(fasta, seqs, out):
    ''' Read bnry file and calculate the category of each cluster '''

    with open(fasta, 'r') as f, open(out, 'w') as o:
        for name, seq in read_fasta(f):
            name = name.split(' ')[0]
            if name in seqs:
                o.write(f"{name}\n{seq}\n")


def get_seqs_in_clstr(clstr, cnumb):
    ''' reads CD-HIT clstr file and retrieves seq names '''

    seqs = {}

    with open(clstr, 'r') as f:

        for l in f:

            if l.startswith('>'):
                cluster = int(l.rstrip().split(' ')[1])

            if not l.startswith('>') and cluster == cnumb:
                gene = l.split(', ')[1].split('...')[0]
                seqs[gene] = ''

    return seqs

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-c', '--cluster_file',
        help='Please specify the CD-HIT Cluster File!',
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-f', '--fasta_file',
        help='Please specify the input fasta file!',
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-n', '--cluster_number',
        help='Please specify the cluster number (ex: 9)!',
        metavar='',
        type=int,
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output fasta file name',
        metavar='',
        type=str,
        )
    args=vars(parser.parse_args())

    clstr = args['cluster_file']
    cnumb = args['cluster_number']
    fasta = args['fasta_file']
    out = args['output_file']

    # Retrieve the list of sequence names for cnumb cluster
    print(f'Retrieving the sequence names for cluster {cnumb}...')
    seqs = get_seqs_in_clstr(clstr, cnumb)

    # Retrieve the fasta sequences for seq names in seqs. Write fasta output.
    print(f'Retrieving fasta sequences for cluster {cnumb} and writing to {out}')
    get_seq_fasta(fasta, seqs, out)


if __name__ == "__main__":
    main()
