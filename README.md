# Pangenome Analysis Pipeline using Prodigal and CD-HIT

This is a new repo. I'm still working on it.

Scripts written for Python version 3.6+


## Step 00: Required tools :: Python 3.6+, Prodigal and CD-HIT.

### Python 3.6+ for running the Python scripts in this repo.

Information for installing and running Python can be found [here](https://www.python.org/). I recommend installing [mini conda](https://docs.conda.io/en/latest/miniconda.html) first and then creating an environment for Python 3.6+ and other tools for the project at hand.

*All Python scripts in this repo were written for Python 3.6+. If you get a syntax error the first time you run a script, please first check your Python version.*

Information for the scripts can be obtained by executing:

```bash
python scriptname.py -h
```

### Prodigal for protein coding gene prediction.
 
Information and installation instructions for Prodigal can be found [here](https://github.com/hyattpd/Prodigal). The publication is [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/).

Prodigal can also be installed with a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```bash
conda create -n prodigal
conda activate prodigal
conda install -c bioconda prodigal
```

### CD-HIT to cluster predicted gene sequence.

Information and installation for CD-HIT can be found [here](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide). The publication can be found [here](https://academic.oup.com/bioinformatics/article/22/13/1658/194225).

CD-HIT can also be installed with a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```bash
conda create -n cdhit
conda activate cdhit
conda install -c bioconda cd-hit
```

## Step 01: Predict protein coding genes with Prodigal.

## Step 02: CD-Hit

## Step 03: Generate binary presence/absence matrix

## Step 04: Build pangenome rarefaction models and curves

Currently, if you've clustered sequences with CD-HIT from multiple genomes, the clstr to binary matrix script will turn the CD-HIT cluster file into a binary matrix with genes as the rows and genomes as the columns. This script needs to be modified at line 68 so the genomeID variable retrieves the unique genomeID from your fasta sequence names.

Once you have a binary matrix, the calculate model plot script will build some nice output like so:

![alt text](Example_Plot.png "Example Pangenome Curve plot.")

## Step 05: Build clustered heatmap based on shared gene content.