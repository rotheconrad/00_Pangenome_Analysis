#!/usr/bin/env python

'''Calculate and Plot a Pangenome Rarefaction Curve from multiple trials

This script reads the B_*org_Rarefaction_Results directory of *.tsv
files from the PGE_0000_Summary folder output by the PGE pipeline.

It models and plots the Pangnome curves from the multiple random trials
(experiments) sampled by the PGE pipeline.

Pangenome modelling based on:
Zhang et. al. 2018, 10.3389/fmicb.2018.00577
Tettelin et. al. 2008, https://doi.org/10.1016/j.mib.2008.09.006

-------------------------------------------
Author :: Roth Conrad & Carlos Ruiz
Email :: rotheconrad@gatech.edu, cruizperez3@gatech.edu
GitHub :: https://github.com/rotheconrad, https://github.com/cruizperez/
Date Created :: January 24th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad & Carlos Ruiz
All rights reserved
-------------------------------------------
'''

import argparse, os
from collections import defaultdict
import pandas as pd
import numpy as np
from lmfit.models import PowerLawModel, ExpressionModel
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_pangenome_curve(
                        dfout, PLM_pan, PLM_spec,
                        omega_core, omega_new,
                        org, prm, trials, out, ymax, ystep
                        ):
    ''' This function builds a plot of the pangenome calculations'''

    print('\n\nBuilding a plot of the data ...')

    # Set the organism name for plot title
    oName = ' '.join(org.split('_'))
    # Set x-axis values
    x = dfout['n'].unique()

    # Calculate mean, median, and quartiles
    dfmean = dfout.groupby('n').mean()
    dfmedian = dfout.groupby('n').median()
    df025 = dfout.groupby('n').quantile(q=0.025)
    df975 = dfout.groupby('n').quantile(q=0.975)

    # Set the colors
    H1 = '#933b41'
    pan = '#933b41'
    core = '#0868ac'
    spec = '#b35806'
    new = '#542788'
    other = '#bdbdbd'
    # alpha value to use for IQRs
    a = 0.2
    # model marker, alpha, and size
    mdl = 'd'
    m = 0.3
    md = 5

    # Build the plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        gridspec_kw={'height_ratios': [2, 1]},
        figsize=(20, 14),
        sharex=True,
        sharey=False
        )

    # plot title, labels, and text
    ax1.text(
        0.5, 1.13, oName,
        fontstyle='italic', fontweight='heavy', fontsize=60,
        color=H1, horizontalalignment='center',
        transform=ax1.transAxes
        )
    ax1.set_title(
        'Pangenome Rarefaction & Core Genes Curve',
        color=H1, fontsize=50, y=1.02
        )
    ax2.set_title(
        'Number of New Genes per Genome / Number of Genes in the Genome',
        color=new, fontsize=32, y=1.02
        )
    ax1.set_ylabel('Number of Gene Clusters', fontsize=28)
    ax2.set_ylabel('New Genes Ratio', fontsize=28)
    ax2.set_xlabel('Number of Genomes', fontsize=28)

    # Emperical data text
    stext = (
        f"Pangenome Size: {dfmean.Pangenome.tolist()[-1]}  |  "
        f"Core Genes: {dfmean.CoreGenome.tolist()[-1]}  |  "
        f"Specific Genes: {dfmean.GenomeSpecific.tolist()[-1]}"
        )
    ax1.text(
        0.5, 0.02, stext,
        fontsize=18, color='#737373',
        horizontalalignment='center', transform=ax1.transAxes
        )
    ttext = (
        f"Number of Genomes: {len(x)+1}  |  "
        f"Number of Permutations: {prm}\n"
        f"Number of Random Trials {trials}"
        )
    ax1.text(
        0.5, 0.90, ttext,
        fontsize=22, color='#737373',
        horizontalalignment='center', transform=ax1.transAxes
        )
    mtext = (
        f"Mean Genes per Genome: {int(dfmean.GenomeLength.mean())}  |  "
        f"Mean New Genes per Genome: {int(dfmean.NewGenes.mean())}  |  "
        f"New Genes at n=100: {int(dfmean.NewGenes.tolist()[-1])}"
        )
    ax2.text(
        0.5, 0.90, mtext,
        fontsize=18, color='#737373',
        horizontalalignment='center', transform=ax2.transAxes
        )

    # Modelled data text
    pangamma = PLM_pan['Gamma']
    specgamma = PLM_spec['Gamma']

    if pangamma < 0: panlabel = 'Closed'
    elif pangamma <= 1: panlabel = 'Open'
    else: print('Error in Pangenome Model Gamma Paramerter')
    if specgamma < 0: speclabel = 'Closed'
    elif specgamma <= 1: speclabel = 'Open'
    else: print('Error in Genome Specific Model Gamma Paramerter')

    modeltext = (
        f'Pangenome Model \u03B3 = {pangamma:.2f}, {panlabel}  |  '
        f'Genome-Specific Model \u03B3 = {specgamma:.2f}, {speclabel}  |  '
        f'Core Gene Model \u03A9 = {omega_core:.2f}'
        )

    ax1.text(
        0.5, -0.05, modeltext,
        fontsize=18, color='#000000',
        horizontalalignment='center', transform=ax1.transAxes
        )
    rtext = (
        f'Ratio at n=100: {dfmean.NewGeneRatio.tolist()[-1]:.2f}%\n'
        f'New Gene Ratio Model \u03A9 = {omega_new:.2f}%'
        )
    ax2.text(
        len(x)-2, dfmean.NewGeneRatio.mean() + 5, rtext,
        fontsize=18, color=new, horizontalalignment='right'
        )


    # Plot Pangenome Medians, Means and IQR
    ax1.plot(x, dfmedian.Pangenome, color=pan, linestyle='--', lw=1)
    ax1.plot(x, dfmean.Pangenome, color=pan, linestyle=':', lw=2)
    ax1.fill_between(x, df025.Pangenome, df975.Pangenome, color=pan, alpha=a)
    ax1.plot(
        x, dfmean.Pangenome_PLM,
        color=pan, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # Plot Genome Specific Medians, Means and IQR
    ax1.plot(x, dfmedian.GenomeSpecific, color=spec, linestyle='--', lw=1)
    ax1.plot(x, dfmean.GenomeSpecific, color=spec, linestyle=':', lw=2)
    ax1.fill_between(
                    x, df025.GenomeSpecific, df975.GenomeSpecific,
                    color=spec, alpha=a
                    )
    ax1.plot(x, dfmean.GenomeSpecific_PLM,
        color=spec, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # Plot Core Medians, Means and IQR
    ax1.plot(x, dfmedian.CoreGenome, color=core, linestyle='--', lw=1)
    ax1.plot(x, dfmean.CoreGenome, color=core, linestyle=':', lw=2)
    ax1.fill_between(x, df025.CoreGenome, df975.CoreGenome, color=core, alpha=a)
    ax1.plot(x, dfmean.CoreGenome_EDM,
        color=core, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # Build ax1 Plot Legend
    l_pan = Line2D(
        [0],[0], color='w', label='Pangenome',
        markerfacecolor=pan, marker='o', markersize=18, alpha=m
        )
    l_core = Line2D(
        [0],[0], color='w', label='Core Genes',
        markerfacecolor=core, marker='o', markersize=18, alpha=m
        )
    l_specific = Line2D(
        [0],[0], color='w', label='Genome-Specific Genes',
        markerfacecolor=spec, marker='o', markersize=18, alpha=m
        )
    l_IQR = Line2D(
        [0],[0], color='w', label='95% IQR',
        markerfacecolor=other, marker='s', markersize=20
        )
    l_mean = Line2D(
        [0],[0], color=other, linestyle=':', lw=4, label='Mean'
        )
    l_median = Line2D(
        [0],[0], color=other, linestyle='--', lw=4, label='Median'
        )
    l_model = Line2D(
        [0],[0], color='w', label='Model Fit',
        markerfacecolor=other, marker=mdl, markersize=20
        )

    legend_elements = [
                        l_pan,
                        l_core,
                        l_specific,
                        l_IQR,
                        l_mean,
                        l_median,
                        l_model
                        ]

    ax1.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=18,
        fancybox=True,
        framealpha=0.0,
        frameon=False
        )

    # Plot Gene Ratio Plot
    ax2.plot(x, dfmedian.NewGeneRatio, color=new, linestyle='--', lw=1)
    ax2.plot(x, dfmean.NewGeneRatio, color=new, linestyle=':', lw=2)
    ax2.fill_between(
                x, df025.NewGeneRatio, df975.NewGeneRatio,
                color=new, alpha=a
                )
    ax2.plot(x, dfmean.NewGeneRatio_EDM,
        color=new, marker=mdl, linestyle='', markersize=md, alpha=m
        )

    # set the axis parameters / style
    ax1.yaxis.grid(which="both", color='#d9d9d9', linestyle='--', linewidth=1)
    ax1.set_yticks(range(0, ymax+1, ystep))
    ax1.minorticks_on()
    ax1.set_xlabel('')
    ax1.tick_params(labelsize=22)
    for spine in ax1.spines.values(): spine.set_linewidth(2)
    ax1.set_axisbelow(True)

    ax2.yaxis.grid(which="both", color='#d9d9d9', linestyle='--', linewidth=1)
    ax2.minorticks_on()
    ax2.tick_params(labelsize=22)
    ax2.set_xticks(range(0, len(x)+2, 10))
    ax2.set_xlim(-1, len(x)+2)
    for spine in ax2.spines.values(): spine.set_linewidth(2)
    ax2.set_axisbelow(True)

    plt.subplots_adjust(
        left = 0.09,
        right = 0.98,
        bottom = 0.07,
        top = 0.87,
        hspace = 0.2
        )

    plt.savefig(f'{out}_pangenome_curves.png')
    plt.close()


def model_decay_curve(dfout, Column):
    ''' This function models the decay curve for core, and new genes'''

    # Model the Core gene curve using an Exponential Decay Function:
    # Fc = Kc*exp(-N/τc) + Ω
    print(f'\nFitting Exponential Decay function to {Column}')
    print('Using Exponential Decay Function: K*exp-(N/\u03C4) + \u03A9')
    # Initialize model
    Custom_Exponential = ExpressionModel(
                                    'A * exp(-x/tau) + omega',
                                    independent_vars=['x']
                                    )
    # Initialize custom parameters
    Expression_Params = Custom_Exponential.make_params()
    # add params with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
    Expression_Params.add_many(
                            ('A', 5, True, 0, None, None, 0.1),
                            ('tau', 5, True, 0, None, None, 0.1),
                            ('omega', 5, True, 0, None, None, 0.1)
                            )
    EDM_Fit = Custom_Exponential.fit(
                                    dfout[Column],
                                    Expression_Params,
                                    x=dfout['n']
                                    )
    dfout[f'{Column}_EDM'] = EDM_Fit.best_fit
    omega = EDM_Fit.best_values['omega']

    return dfout, omega


def model_pangenome_curve(dfout, Column):
    ''' This function models the panganome, core, and new gene curves'''

    # Initialize dictionary to store model results
    params = {} # PowerLawModel for Pangenome curve
    # Model the Pangenome curve using a Powerlaw function: Ps = κn^γ
    print(f'\nFitting PowerLaw function to {Column}')
    print('Using Power Law Function K*N^\u03B3 ...')
    PLM = PowerLawModel() # Initialize the model
    # Guess starting values from the data
    PLM_Guess_Start = PLM.guess(
                            dfout[Column],
                            x=dfout['n'],
                            )
    # Fit model to the data (optimize parameters)
    PLM_Fit = PLM.fit(
                    dfout[Column],
                    PLM_Guess_Start,
                    x=dfout['n']
                    )
    # Add modeled data to dfout data frame
    dfout[f'{Column}_PLM'] = PLM_Fit.best_fit

    # Store K and Gamma parameters 
    params['K'] = float(PLM_Fit.best_values['amplitude'])
    params['Gamma'] = float(PLM_Fit.best_values['exponent'])

    return dfout, params

def collect_results(fdir, out):
    ''' reads directory of tsv files, builds combined dataframe. '''

    # Define columns for output data frame
    colout = [
            'Trial', 'Permutation', 'n', 'Pangenome', 'CoreGenome',
            'GenomeSpecific', 'NewGenes', 'NewGeneRatio', 'GenomeLength'
            ]
    # Initialize list to store new rows
    data = []
    # Create list of files from tsv directory
    if fdir[-1] == '/': fdir = fdir[:-1]
    f_list = [f for f in os.listdir(fdir) if os.path.isfile(f'{fdir}/{f}')]
    # define number of trials as number of files in fdir
    trials = len(f_list)

    # Read through each file and transform the tsv into extended dictionary
    for trial, file in enumerate(f_list):
        # Track trial number = to file number
        trial_number = trial + 1
        print(f'Collecting results data from file 1 of {trial_number} ...')
        # Read results tsv and reconstruct dataframe
        with open(f'{fdir}/{file}', 'r') as f:
            skip_header = f.readline()
            skip_first_genome = f.readline()
            n = 2
            for line in f:
                # initialize empty dict for each line
                d = {}
                X = line.rstrip().split('\t')
                d['Pangenome'] = [int(x) for x in X[6][1:-1].split(', ')]
                d['CoreGenome'] = [int(x) for x in X[2][1:-1].split(', ')]
                d['GenomeSpecific'] = [int(x) for x in X[8][1:-1].split(', ')]
                d['NewGenes'] = [int(x) for x in X[4][1:-1].split(', ')]
                d['NewGeneRatio'] = [float(x) for x in X[5][1:-1].split(', ')]
                d['GenomeLength'] = [int(x) for x in X[3][1:-1].split(', ')]
                steps = len(d['Pangenome'])

                for i in range(steps):
                    newrow = [
                            trial_number,
                            i+1,
                            n,
                            d['Pangenome'][i],
                            d['CoreGenome'][i],
                            d['GenomeSpecific'][i],
                            d['NewGenes'][i],
                            d['NewGeneRatio'][i],
                            d['GenomeLength'][i]
                            ]
                    data.append(newrow)

                n += 1

    dfout = pd.DataFrame(data=data, columns=colout)
    # define permutations as length of one of the lists
    prms = len(X[6][1:-1].split(', '))

    return(dfout, prms, trials)

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--results_directory',
        help='Please specify the binary.tsv input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-n', '--organism_name',
        help='Organism_name_experiment ex: Escherichia_coli_0042',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-y', '--yaxis_max',
        help='Set range of y-axis pangenome size (eg: 20000)',
        metavar='',
        type=int,
        required=True
        )
    parser.add_argument(
        '-s', '--yaxis_step',
        help='Set y-axis step increment (eg: 2000)',
        metavar='',
        type=int,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='What do you want to name the output files?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    print('\n\nRunning Script ...\n\n')
    # Collect the result data from the tsv files
    dfout, prms, trials = collect_results(
                            args['results_directory'],
                            args['output_prefix']
                            )

    # Fit a Power Law Model to the Pangenome Curve
    dfout, PLM_pan = model_pangenome_curve(dfout, 'Pangenome')
    # Fit a Power Law Model to the Genome Specific Genome
    dfout, PLM_spec = model_pangenome_curve(dfout, 'GenomeSpecific')
    # Fit an Exponential Decay Model to the Core Genome
    dfout, omega_core = model_decay_curve(dfout, 'CoreGenome')
    # Fit an Exponential Decay Model to the New Genes per genome count
    dfout, omega_new = model_decay_curve(dfout, 'NewGeneRatio')
    # Plot the results
    plot_pangenome_curve(
                        dfout,
                        PLM_pan,
                        PLM_spec,
                        omega_core,
                        omega_new,
                        args['organism_name'],
                        prms,
                        trials,
                        args['output_prefix'],
                        args['yaxis_max'],
                        args['yaxis_step']
                        )

    print(
        '\n\nCongratulations!! '
        'The script seems to have finished successfully.\n\n'
        )

if __name__ == "__main__":
    main()