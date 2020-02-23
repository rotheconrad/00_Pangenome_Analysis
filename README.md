# Pangenome Analysis Pipeline using Prodigal and CD-HIT

This is a new repo. I'm still working on it.

Scripts written for Python version 3.6+

Information for the scripts can be obtained by executing:

```bash
python scriptname.py -h
```

Currently, if you've clustered sequences with CD-HIT from multiple genomes, the clstr to binary matrix script will turn the CD-HIT cluster file into a binary matrix with genes as the rows and genomes as the columns. This script needs to be modified at line 68 so the genomeID variable retrieves the unique genomeID from your fasta sequence names.

Once you have a binary matrix, the calculate model plot script will build some nice output like so:

![alt text](Example_Plot.png "Example Pangenome Curve plot.")