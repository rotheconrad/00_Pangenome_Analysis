#!/usr/bin/env python

''' Append Unique Genome or MAG ID to beginning of contig names

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

import argparse, subprocess

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

def Fasta_rename_sequences(infile, uniqueid):

    outfile = infile + '.rename'

    with open(infile, 'r') as f, open(outfile, 'w') as o:

        for name, seq in read_fasta(f):
            pname = f'>{uniqueid}_{name[1:]}'
            newName = f'{pname}\n{seq}\n'
            o.write(newName)

    _ = subprocess.run(['mv', outfile, infile])


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify an input fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-u', '--unique_id',
        help='Please specify the unique ID to append!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')
    infile = args['input_file']
    uniqueid = args['unique_id']
    Fasta_rename_sequences(infile, uniqueid)
    
if __name__ == "__main__":
    main()

