#!/usr/bin/env python

''' Filters fasta file for minimum sequence length

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

def Fasta_rename_sequences(infile, len_filter):

    outfile = infile + '.lenfilter'

    data = []

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        i = 0
        c = 0
        for name, seq in read_fasta(f):
            c += 1
            if len(seq) > len_filter:
                o.write(f'{name}\n{seq}\n')
                i += 1

            elif len(seq) <= len_filter:
                data.append(len(seq))

    _ = subprocess.run(['mv', outfile, infile])

    print(f'\n\n#### {infile} ################################################')
    print(f'Kept {i} of {c} predicted gene sequences')
    print(f'Sequences below {len_filter} nucleotide filter cutoff: {c-i}')
    print(f'Shortest gene length: {min(data)}')
    print(f'Total genes filtered: {len(data)}')
    print(f'##############################################################\n\n')

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
        '-m', '--minimum_length',
        help='Please specify the unique ID to append!',
        metavar='',
        type=int,
        required=True
        )
    args=vars(parser.parse_args())

    infile = args['input_file']
    len_filter = args['minimum_length']
    Fasta_rename_sequences(infile, len_filter)
    
if __name__ == "__main__":
    main()

