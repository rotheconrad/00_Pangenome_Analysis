#!/usr/local/pacerepov1/python/2.7/bin/python

## USAGE :: python scriptname.py file.fasta prefix
## Reads fasta file and replaces sequence names with prefix_# counting from 1 to the end.
## Writes file.ref matching original names to new names.

import sys, subprocess

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

def Fasta_rename_sequences(infile, minlen):

#    X = infile.split('/')
#    prefix = X[1].split('.')[0]
#    outfile = X[0] + prefix + '.rename'
    prefix = infile.split('.')[0]
    outfile = infile + '.rename'

    with open(infile, 'r') as f, open(outfile, 'w') as o:

        for name, seq in read_fasta(f):
            pname = name[1:].split(' ')[0]
            length = int(pname.split('_')[3])
            if length > minlen:
                newName = f'>{prefix}_{pname}\n{seq}\n'
                o.write(newName)


    _ = subprocess.run(['mv', outfile, infile])


def main():
    infile = sys.argv[1]
    minlen = int(sys.argv[2])
    Fasta_rename_sequences(infile, minlen)
    
if __name__ == "__main__":
    main()

