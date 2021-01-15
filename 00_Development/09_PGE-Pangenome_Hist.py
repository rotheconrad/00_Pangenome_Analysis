#!~/.conda/envs/rothon/bin/python

## USAGE :: python scriptname.py -cd Comparison_Directory -n number_of_genomes_compared -op output_prefix
## Reads through directory of pangenome comparison files and builds two summary plots.
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Saturday, June 22, 2019
## Date Updated :: N/A first version

def plot_pangenome_summary(df, op):
    ''' Takes a data frame with columns as species names and rows and pangenome statistics. plots bar plots'''

    grid_color, xgrid_color ='#d9d9d9', '#bdbdbd'
    cpan1, cpan2 = '#933b41', '#ddb7b1'
    cnew1, cnew2 = '#c9d7e9', '#a6bddb'
    cspecific, crare, ccommon, ccore = '#fdb863', '#b2abd2', '#41ab5d', '#4eb3d3'
    mean_line = "#000000"
    median_line = "#000000"
    alpha = 0.6
    vline_alpha = 0.25
    bins = 10

    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16)) = plt.subplots(4, 4, figsize=(20, 14), sharex=False, sharey=False, gridspec_kw={'height_ratios': [1.6, 1, 1, 1]})

    sepline = '-- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- --'
    panfrac ='-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  Number of Genes as Fraction of Pangenome Size  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    sizefrac = '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- --  Number of Genes as Fraction of Genome Size  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    plt.text(1.09, -0.11, sepline, color='#737373', fontsize=16, fontweight='heavy', horizontalalignment='center', transform=ax2.transAxes)
    plt.text(1.08, 1.1, panfrac, color='#737373', fontsize=16, fontweight='heavy', horizontalalignment='center', transform=ax10.transAxes)
    plt.text(1.08, 1.1,  sizefrac, color='#737373', fontsize=16, fontweight='heavy', horizontalalignment='center', transform=ax14.transAxes)


    ## Row 1 - Pangenome and New Genes
    ax1.hist(df.Pangenome, color=cpan1, alpha=alpha, bins=bins)
    ax1.set_title('Pangenome Size', fontweight='heavy', color=cpan1, fontsize=20, y=1.02)
    ax1.set_ylabel('Genes', fontweight='bold', fontsize=14)
    ax1.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    ax1.minorticks_on()
    ax1.set_axisbelow(True)

    ax2.hist(df.Pan_norm, color=cpan2, alpha=alpha, bins=bins)
    ax2.set_title('Pangenome / Genome Size', fontweight='heavy', color=cpan2, fontsize=20, y=1.02)
    ax2.set_ylabel('Genes / Genome Size', fontweight='bold', fontsize=14)
    ax2.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax2.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax2.minorticks_on()
    ax2.set_axisbelow(True)

    ax3.hist(df.New_Genes, color=cnew1, alpha=alpha, bins=bins)
    ax3.set_title('New Genes per Genome', fontweight='heavy', color=cnew1, fontsize=20, y=1.02)
    ax3.set_ylabel('Genes', fontweight='bold', fontsize=14)
    ax3.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax3.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax3.minorticks_on()
    ax3.set_axisbelow(True)

    ax4.hist(df.New_norm, color=cnew2, alpha=alpha, bins=bins)
    ax4.set_title('New Genes / Genome Size', fontweight='heavy', color=cnew2, fontsize=20, y=1.02)
    ax4.set_ylabel('Genes / Genome Size', fontweight='bold', fontsize=14)
    ax4.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax4.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax4.minorticks_on()
    ax4.set_axisbelow(True)

    # Row 2 - Gene Counts
    ax5.hist(df.Specific, color=cspecific, alpha=alpha, bins=bins)
    ax5.set_title('Specific Genes', fontweight='heavy', color=cspecific, fontsize=24)
    ax5.set_ylabel('Gene Counts', fontweight='bold', fontsize=14)
    ax5.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax5.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax5.minorticks_on()
    ax5.set_axisbelow(True)
    ax5ytick = ax5.get_yticks()
    ax5.set_yticks(ax5ytick)


    ax6.hist(df.Rare, color=crare, alpha=alpha, bins=bins)
    ax6.set_title('Rare Genes', fontweight='heavy', color=crare, fontsize=24)
    ax6.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax6.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax6.minorticks_on()
    ax6.set_axisbelow(True)
    ax6.set_yticks(ax5ytick)

    ax7.hist(df.Common, color=ccommon, alpha=alpha, bins=bins)
    ax7.set_title('Common Genes', fontweight='heavy', color=ccommon, fontsize=24)
    ax7.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax7.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax7.minorticks_on()
    ax7.set_axisbelow(True)
    ax7.set_yticks(ax5ytick)

    ax8.hist(df.Core, color=ccore, alpha=alpha, bins=bins)
    ax8.set_title('Core Genes', fontweight='heavy', color=ccore, fontsize=24)
    ax8.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax8.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax8.minorticks_on()
    ax8.set_axisbelow(True)
    ax8.set_yticks(ax5ytick)

    # Row 3 - Genes / Pangenome
    df.Specific_frac.plot.hist(color=cspecific, ax=ax9, alpha=alpha, bins=bins)
    ax9.set_ylabel('Genes / Pangenome', fontweight='bold', fontsize=14)
    ax9.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax9.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax9.minorticks_on()
    ax9.set_axisbelow(True)
    ax9ytick = ax9.get_yticks()
    ax9.set_yticks(ax9ytick)

    df.Rare_frac.plot.hist(color=crare, ax=ax10, alpha=alpha, bins=bins)
    ax10.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax10.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax10.minorticks_on()
    ax10.set_axisbelow(True)
    ax10.set_yticks(ax9ytick)

    df.Common_frac.plot.hist(color=ccommon, ax=ax11, alpha=alpha, bins=bins)
    ax11.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax11.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax11.minorticks_on()
    ax11.set_axisbelow(True)
    ax11.set_yticks(ax9ytick)

    df.Core_frac.plot.hist(color=ccore, ax=ax12, alpha=alpha, bins=bins)
    ax12.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax12.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax12.minorticks_on()
    ax12.set_axisbelow(True)
    ax12.set_yticks(ax9ytick)

    # Row 4 - Genes / Genome Size
    df.Specific_norm.plot.hist(color=cspecific, ax=ax13, alpha=alpha, bins=bins)
    ax13.set_ylabel('Genes / Genome Size', fontweight='bold', fontsize=14)
    ax13.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax13.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax13.minorticks_on()
    ax13.set_axisbelow(True)
    ax13ytick = ax13.get_yticks()
    ax13.set_yticks(ax13ytick)

    df.Rare_norm.plot.hist(color=crare, ax=ax14, alpha=alpha, bins=bins)
    ax14.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax14.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax14.minorticks_on()
    ax14.set_axisbelow(True)
    ax14.set_yticks(ax13ytick)

    df.Common_norm.plot.hist(color=ccommon, ax=ax15, alpha=alpha, bins=bins)
    ax15.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax15.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax15.minorticks_on()
    ax15.set_axisbelow(True)
    ax15.set_yticks(ax13ytick)

    df.Core_norm.plot.hist(color=ccore, ax=ax16, alpha=alpha, bins=bins)
    ax16.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax16.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax16.minorticks_on()
    ax16.set_axisbelow(True)
    ax16.set_yticks(ax13ytick)

    # Mean / Median Lines
    ax1.axvline(x=df.Pangenome.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax1.axvline(x=df.Pangenome.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax2.axvline(x=df.Pan_norm.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax2.axvline(x=df.Pan_norm.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax3.axvline(x=df.New_Genes.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax3.axvline(x=df.New_Genes.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax3.axvline(x=df.New_Genes.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax3.axvline(x=df.New_Genes.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax4.axvline(x=df.New_norm.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha, label="Mean")
    ax4.axvline(x=df.New_norm.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha, label="Median")
    ax4.legend(loc='upper right', fontsize=18, frameon=False)

    ax5.axvline(x=df.Specific.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax5.axvline(x=df.Specific.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax6.axvline(x=df.Rare.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax6.axvline(x=df.Rare.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax7.axvline(x=df.Common.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax7.axvline(x=df.Common.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax8.axvline(x=df.Core.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax8.axvline(x=df.Core.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax9.axvline(x=df.Specific_frac.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax9.axvline(x=df.Specific_frac.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax10.axvline(x=df.Rare_frac.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax10.axvline(x=df.Rare_frac.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax11.axvline(x=df.Common_frac.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax11.axvline(x=df.Common_frac.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax12.axvline(x=df.Core_frac.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax12.axvline(x=df.Core_frac.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax13.axvline(x=df.Specific_norm.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax13.axvline(x=df.Specific_norm.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax14.axvline(x=df.Rare_norm.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax14.axvline(x=df.Rare_norm.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax15.axvline(x=df.Common_norm.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax15.axvline(x=df.Common_norm.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    ax16.axvline(x=df.Core_norm.mean(), ymin=0, ymax=1, color=mean_line, linewidth=2, linestyle='--', alpha=vline_alpha)
    ax16.axvline(x=df.Core_norm.median(), ymin=0, ymax=1, color=median_line, linewidth=2, linestyle=':', alpha=vline_alpha)

    # Adjust, save, and close
    plt.subplots_adjust(bottom=0.06, hspace=0.3, top=0.94, left=0.09, right=0.96)
    plt.savefig(f'{op}_pangenome_hist.png')
    plt.close()


def parse_pangenome_summary_file(cd, file, pangenome_summary):
    ''' parses a pangenome summary file, populates the dicts and returns them '''

    with open(f'{cd}{file}', 'r') as f:
        # Each pangenome summary file should only be one tab separated line
        X = f.readline().rstrip().split('\t')
        # Select the mean number of genes and pangenome for each organism
        p = float(X[2])
        m = float(X[7])
        # Make a list of the selected variables:
        # X[0] is Organism name. X[1] is Experiment number. X[2] is Pangenome count, X[3] is core genome, X[4] common,
        # X[5] rare, X[6] specific, X[7] mean Genome Size (genes per genome), X[8] New genes (mean new genes per genome)
        name = '_'.join(X[0:2])
        pangenome_summary[name] = [
                            float(X[2]), float(X[2])/m, float(X[8]), float(X[8])/m,
                            float(X[6])/p, float(X[5])/p, float(X[4])/p, float(X[3])/p,
                            float(X[6])/m, float(X[5])/m, float(X[4])/m, float(X[3])/m,
                            float(X[6]), float(X[5]), float(X[4]), float(X[3])
                            ]

    return pangenome_summary


def compare_pangenomes(cd, op):
    ''' Reads the pangenome comparison files in cd, builds two dictionaries, and passes them to functions to build plots '''

    # Initialize Dictionary: order matches list l from parse_pangenome_summary_file() function.
    pangenome_summary = {'Order': ['Pangenome', 'Pan_norm', 'New_Genes', 'New_norm',
                                   'Specific_frac', 'Rare_frac', 'Common_frac', 'Core_frac',
                                   'Specific_norm', 'Rare_norm', 'Common_norm', 'Core_norm',
                                   'Specific', 'Rare', 'Common', 'Core']}

    file_list = os.listdir(cd)

    # Read through files and populate dictionaries
    for file in file_list: pangenome_summary = parse_pangenome_summary_file(cd, file, pangenome_summary)

    # Pangenome Summary Plot
    df = pd.DataFrame(pangenome_summary).set_index(['Order']).T
    df = df.sort_values(by=['Pan_norm'], ascending=False)

    print(df.describe())

    plot_pangenome_summary(df, op)
    #print(df.index)

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(description='Reads through directory of pangenome comparison files and builds two summary plots. Designed to compary 10 or fewer organisms.')
    parser.add_argument('-cd', '--comparison_directory', help='Please specify the directory of pangenome comparison files!', required=True)
    parser.add_argument('-op', '--out_file', help='What do you want to name the output files?', required=True)
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')
    compare_pangenomes(args['comparison_directory'], args['out_file'])

if __name__ == "__main__":
    import argparse, os, matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    #import matplotlib.patches as mpatches
    import matplotlib.lines as lines
    import pandas as pd
    main()
