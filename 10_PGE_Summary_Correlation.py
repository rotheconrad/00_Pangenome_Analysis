#!/usr/bin/env python

'''PGE Summary - ANI vs Shared Fragment Correlation.

This script reads a directory of PGE Correlations summary tsv files
and builds a scatter plot of ANI values on the x-axis and the ratio
of shared fragments over total fragments on the y-axis. It plots in
gray the results of the total PGE results and in black the result of
the single specified representative experiment.

This tool takes the following input parameters:

    * csd - correlation summary directory (str)
    * org - organism name for the plot title (str)
    * rep - representative experiment number (int, ex: 0042)
    * op - An output file prefix (str)

This script returns the following files:

    * .png file of plot image with {op}_correlation_plot.png

This script requires the following packages:

    * argparse
    * os
    * numpy
    * pandas
    * matplotlib

This file can also be imported as a module and contains the following 
functions:

    * ANI_vs_shared_plot - plots the data. writes the png.
    * gather_data - reads files from input directory. Builds DataFrame
    * gather_stats - computes summary stats and correlation on the data
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, June 27th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def gather_stats(df):
    """Computes correlation, mean, and median on df columns xs and ys

    Parameters
    ----------
    df : DataFrame
        DataFrame with columns xs: ANI values, ys: ratio shared fragments

    Returns
    -------
    df_stats : dict
        dictionary with keys scorr, pcorr, kcorr, ani_mean, ani_median,
        frag_mean, frag_median
    """

    # Compute Pearson Correlation Coefficient
    pcorr = round(df['xs'].corr(df['ys'], method='pearson'), 3)
    # Compute Spearman's Correlation
    scorr = round(df['xs'].corr(df['ys'], method='spearman'), 3)
    # Compute Kendall's Ta3
    kcorr = round(df['xs'].corr(df['ys'], method='kendall'), 3)
    # Compute ANI mean and median
    ani_mean = np.mean(df['xs'])
    ani_median = np.median(df['xs'])
    # Compute shared fragment mean and median
    frag_mean = np.mean(df['ys'])
    frag_median = np.median(df['ys'])
    # Compile dictionairy
    df_stats = {
        'pcorr': pcorr,
        'scorr': scorr,
        'kcorr': kcorr,
        'ani_mean': ani_mean,
        'ani_median': ani_median,
        'frag_mean': frag_mean,
        'frag_median': frag_median
        }

    print(f"ANI mean: {ani_mean}", f"ANI median: {ani_median}")
    print(f"Frag mean: {frag_mean}", f"Frag median: {frag_median}")

    return df_stats


def ANI_correlation_plot(df_full, df_rep, org, op):
    """Takes the data and builds the plot

    Parameters
    ----------
    df_full : DataFrame
        Dataframe with data for all the experients run
    df_rep : DataFrame
        Dataframe with data for the selected representative experiment.
    org : str
        The name of the organism for the plot title
    op : str
        The prefix for the output plot.png

    Returns
    -------
    No return. Saves the plot to {op}_correlation_plot.png in the
    working directory unless another directory was specified as
    part of the op.
    """

    org_name = ' '.join(org.split('_'))

    # Gather Stats
    df_full_stats = gather_stats(df_full)
    df_rep_stats = gather_stats(df_rep)
    stats_line = (
        "Pearson Correlation Coefficient:\n"
        "Spearman's Correlation:\n"
        "Kendall's Tau:\n"
        )
    full_stats_line = (
        f"{df_full_stats['pcorr']},\n"
        f"{df_full_stats['scorr']},\n"
        f"{df_full_stats['kcorr']},"
        )
    rep_stats_line = (
        f"{df_rep_stats['pcorr']}\n"
        f"{df_rep_stats['scorr']}\n"
        f"{df_rep_stats['kcorr']}"
        )

    # Set Colors
    full_color = "#bdbdbd"
    full_line_color = "#000000"
    rep_color = "#000000"

    # Build the plot
    fig, ax = plt.subplots(figsize=(20, 10))

    # plot title, labels, and text
    ax.text(
        0.5, 1.13, org_name,
        fontstyle='italic', fontweight='heavy', fontsize=60, color=rep_color,
        horizontalalignment='center', transform=ax.transAxes
        )
    ax.set_title(
        'ANI Value vs Shared Sequence Fragment Ratio',
        fontsize=50, y=1.02, color=rep_color
        )
    ax.text(
        0.24, 0.99, stats_line,
        fontsize=18, color='#252525',
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes
        )
    ax.text(
        0.245, 0.99, f'{full_stats_line}',
        fontsize=18, fontweight='bold', color=full_color,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes
        )
    ax.text(
        0.30, 0.99, f'{rep_stats_line}',
        fontsize=18, fontweight='bold', color=rep_color,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes
        )
    ax.set_xlabel(
        'Average Nucleotide Identity (fastANI)',
        fontsize=30, fontweight='bold', y=-0.02
        )
    ax.set_ylabel(
        'Shared / Total Fragments',
        fontsize=30, fontweight='bold', x=-0.02
        )

    # set the axis parameters / style
    ax.minorticks_on()
    ax.set_xticks(np.arange(96.9, 100.1, 0.1), minor=True)
    ax.set_xticks(np.arange(97.0, 100.1, 0.5))
    ax.set_yticks(np.arange(0.46, 1.02, 0.02), minor=True)
    ax.set_yticks(np.arange(0.5, 1.1, 0.1))
    ax.set_xlim(left=96.9, right=100.1)
    ax.set_ylim(bottom=0.46, top=1.02)
    ax.tick_params(axis='both', labelsize=18)
    ax.tick_params(
        axis='x', which='major', direction='in', color='k',
        width=6, length=12, bottom=True, zorder=3
        )

    """
    # For Pcoccus
    ax.set_xticks(np.arange(86.9, 100.1, 0.1), minor=True)
    ax.set_xticks(np.arange(87.0, 100.1, 0.5))
    ax.set_yticks(np.arange(0.0, 1.02, 0.02), minor=True)
    ax.set_yticks(np.arange(0.0, 1.1, 0.1))

    ax.set_xlim(left=87.9, right=100.1)
    ax.set_ylim(bottom=0.0, top=1.02)
    ax.tick_params(axis='both', labelsize=18)
    ax.tick_params(
        axis='x', which='major', direction='in', color='k',
        width=6, length=12, bottom=True, zorder=3
        )
    """

    # set grid style
    ax.yaxis.grid(which="minor", color='#d9d9d9', linestyle='--', linewidth=1)
    ax.xaxis.grid(which="minor", color='#f0f0f0', linestyle='-', linewidth=1)
    ax.yaxis.grid(which="major", color='#d9d9d9', linestyle='--', linewidth=1.5)
    ax.xaxis.grid(which="major", color='#f0f0f0', linestyle='-', linewidth=2)
    ax.set_axisbelow(True)

    # plot the data for df_full
    ax.scatter(
        df_full['xs'], df_full['ys'],
        marker='o', s=30, color=full_color, alpha=0.15
        )
    full_mn = ax.axvline(
            x=df_full_stats['ani_mean'], ymin=0, ymax=1,
            color=full_line_color, linewidth=3, linestyle='--',
            label='Mean'
            )
    full_md = ax.axvline(
            x=df_full_stats['ani_median'], ymin=0, ymax=1,
            color=full_line_color, linewidth=3, linestyle=':',
            label='Median'
            )
    _ = plt.axhline(
        y=df_full_stats['frag_mean'], xmin=0, xmax=1,
        color=full_line_color, linewidth=3, linestyle='--',
        label='Mean'
        )
    _ = plt.axhline(
        y=df_full_stats['frag_median'], xmin=0, xmax=1,
        color=full_line_color, linewidth=3, linestyle=':',
        label='Median'
        )

    # plot the data for df_rep
    ax.scatter(
        df_rep['xs'], df_rep['ys'],
        marker='o', s=30, color=rep_color, alpha=0.15
        )

    # plot the legend
    plt.legend(
        handles=[full_mn, full_md],
        loc='lower right',
        fontsize=24,
        frameon=False,
        ncol=2
        )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(f'{op}_correlation_plot.png')
    plt.close()    


def gather_data(csd, rep):
    """Gets the data and builds a DataFrame

    Parameters
    ----------
    csd : str
        The directory of PGE Correlation Summary tsv files
    rep : str
        The experiment number to use as the representative data

    Returns
    -------
    df_full : DataFrame
        DataFrame with two columns 'xs' and 'ys' containing all data
    df_rep : DataFrame
        DataFrame with tow columns 'xs' and 'ys' containing only data
        from the specified representative experiment
    """

    data_full = {'xs': [], 'ys': []}
    data_rep = {'xs': [], 'ys': []}

    file_list = os.listdir(csd)

    for file in file_list:
        
        experiment_number = file.split('_')[1]

        try:
            with open(f'{csd}{file}', 'r') as f:
                # We do not need the first line here
                _ = f.readline()
                # Second line contains ANI values for x-axis
                xs = f.readline().rstrip().split('\t')
                # Third line contains shared fragment ratio for y-axis
                ys = f.readline().rstrip().split('\t')

                data_full['xs'].extend([float(x) for x in xs])
                data_full['ys'].extend([float(y) for y in ys])

                if experiment_number == rep:
                    data_rep['xs'].extend([float(x) for x in xs])
                    data_rep['ys'].extend([float(y) for y in ys])
        except:
            print(file)

    df_full = pd.DataFrame(data_full)
    df_full = df_full[df_full['xs'] != 100.0]
    df_rep = pd.DataFrame(data_rep)
    df_rep = df_rep[df_rep['xs'] != 100.0]

    print(df_full.describe())
    print(df_rep.describe())

    return df_full, df_rep


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-csd', '--correlation_summary_directory',
        help='Please specify the PGE correlation summary directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-org', '--organism_name',
        help='Organism name for plot title. ex: Escherichia_coli',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rep', '--representative_experiment_number',
        help='Number of experiment to use as representative. ex: 0042',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-op', '--output_file_prefix',
        help='What do you want to name the ouput file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    df_full, df_rep = gather_data(
        args['correlation_summary_directory'],
        args['representative_experiment_number']
        )
    ANI_correlation_plot(
        df_full,
        df_rep,
        args['organism_name'],
        args['output_file_prefix']
        )

if __name__ == "__main__":
    main()
