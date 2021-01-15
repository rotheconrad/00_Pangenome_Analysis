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
    alpha = 0.6
    capsize = 5

    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16)) = plt.subplots(4, 4, figsize=(20, 14), sharex=True, sharey=False, gridspec_kw={'height_ratios': [1.6, 1, 1, 1]})

    #f.suptitle('Bar Plots Showing Comparisons Between Pangenomes', y=0.95, fontsize=24)
    sepline = '-- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- --'
    panfrac ='-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  Number of Genes as Fraction of Pangenome Size  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    sizefrac = '-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- --  Number of Genes as Fraction of Genome Size  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --'
    plt.text(1.09, -0.1, sepline, color='#737373', fontsize=16, fontweight='heavy', horizontalalignment='center', transform=ax2.transAxes)
    plt.text(1.08, 1.143, panfrac, color='#737373', fontsize=16, fontweight='heavy', horizontalalignment='center', transform=ax10.transAxes)
    plt.text(1.08, 1.143,  sizefrac, color='#737373', fontsize=16, fontweight='heavy', horizontalalignment='center', transform=ax14.transAxes)


    ## Row 1 - Pangenome and New Genes
    df.Pangenome.plot.bar(color=cpan1, width=0.8, ax=ax1, alpha=alpha, yerr=[df.Pan_error_low, df.Pan_error_high], capsize=capsize)
    ax1.set_title('Pangenome Size', fontweight='heavy', color=cpan1, fontsize=20, y=1.02)
    ax1.set_ylabel('Genes', fontweight='bold', fontsize=14)
    ax1.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax1.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax1.minorticks_on()
    ax1.set_axisbelow(True)

    df.Pan_norm.plot.bar(color=cpan2, width=0.8, ax=ax2, alpha=alpha, yerr=[df.Pan_error_low_m, df.Pan_error_high_m], capsize=capsize)
    ax2.set_title('Pangenome / Genome Size', fontweight='heavy', color=cpan2, fontsize=20, y=1.02)
    ax2.set_ylabel('Genes / Genome Size', fontweight='bold', fontsize=14)
    ax2.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax2.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax2.minorticks_on()
    ax2.set_axisbelow(True)

    df.New_Genes.plot.bar(color=cnew1, width=0.8, ax=ax3, alpha=alpha, yerr=[df.New_error_low, df.New_error_high], capsize=capsize)
    ax3.set_title('New Genes per Genome', fontweight='heavy', color=cnew1, fontsize=20, y=1.02)
    ax3.set_ylabel('Genes', fontweight='bold', fontsize=14)
    ax3.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax3.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax3.minorticks_on()
    ax3.set_axisbelow(True)

    df.New_norm.plot.bar(color=cnew2, width=0.8, ax=ax4, alpha=alpha, yerr=[df.New_error_low_m, df.New_error_high_m], capsize=capsize)
    ax4.set_title('New Genes / Genome Size', fontweight='heavy', color=cnew2, fontsize=20, y=1.02)
    ax4.set_ylabel('Genes / Genome Size', fontweight='bold', fontsize=14)
    ax4.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax4.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax4.minorticks_on()
    ax4.set_axisbelow(True)

    # Row 2 - Gene Counts
    df.Specific.plot.bar(color=cspecific, width=0.8, ax=ax5, alpha=alpha, yerr=[df.Spec_error_low, df.Spec_error_high], capsize=capsize)
    ax5.set_title('Specific Genes', fontweight='heavy', color=cspecific, fontsize=24)
    ax5.set_ylabel('Gene Counts', fontweight='bold', fontsize=14)
    ax5.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax5.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax5.minorticks_on()
    ax5.set_axisbelow(True)
    ax5ytick = ax5.get_yticks()
    ax5.set_yticks(ax5ytick)


    df.Rare.plot.bar(color=crare, width=0.8, ax=ax6, alpha=alpha, yerr=[df.Rare_error_low, df.Rare_error_high], capsize=capsize)
    ax6.set_title('Rare Genes', fontweight='heavy', color=crare, fontsize=24)
    ax6.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax6.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax6.minorticks_on()
    ax6.set_axisbelow(True)
    ax6.set_yticks(ax5ytick)

    df.Common.plot.bar(color=ccommon, width=0.8, ax=ax7, alpha=alpha, yerr=[df.Common_error_low, df.Common_error_high], capsize=capsize)
    ax7.set_title('Common Genes', fontweight='heavy', color=ccommon, fontsize=24)
    ax7.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax7.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax7.minorticks_on()
    ax7.set_axisbelow(True)
    ax7.set_yticks(ax5ytick)

    df.Core.plot.bar(color=ccore, width=0.8, ax=ax8, alpha=alpha, yerr=[df.Core_error_low, df.Core_error_high], capsize=capsize)
    ax8.set_title('Core Genes', fontweight='heavy', color=ccore, fontsize=24)
    ax8.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax8.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax8.minorticks_on()
    ax8.set_axisbelow(True)
    ax8.set_yticks(ax5ytick)

    # Row 3 - Genes / Pangenome
    df.Specific_frac.plot.bar(color=cspecific, width=0.8, ax=ax9, alpha=alpha, yerr=[df.Spec_error_low_p, df.Spec_error_high_p], capsize=capsize)
    ax9.set_ylabel('Genes / Pangenome', fontweight='bold', fontsize=14)
    ax9.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax9.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax9.minorticks_on()
    ax9.set_axisbelow(True)
    ax9ytick = ax9.get_yticks()
    ax9.set_yticks(ax9ytick)

    df.Rare_frac.plot.bar(color=crare, width=0.8, ax=ax10, alpha=alpha, yerr=[df.Rare_error_low_p, df.Rare_error_high_p], capsize=capsize)
    ax10.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax10.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax10.minorticks_on()
    ax10.set_axisbelow(True)
    ax10.set_yticks(ax9ytick)

    df.Common_frac.plot.bar(color=ccommon, width=0.8, ax=ax11, alpha=alpha, yerr=[df.Common_error_low_p, df.Common_error_high_p], capsize=capsize)
    ax11.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax11.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax11.minorticks_on()
    ax11.set_axisbelow(True)
    ax11.set_yticks(ax9ytick)

    df.Core_frac.plot.bar(color=ccore, width=0.8, ax=ax12, alpha=alpha, yerr=[df.Core_error_low_p, df.Core_error_high_p], capsize=capsize)
    ax12.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax12.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax12.minorticks_on()
    ax12.set_axisbelow(True)
    ax12.set_yticks(ax9ytick)

    # Row 4 - Genes / Genome Size
    df.Specific_norm.plot.bar(color=cspecific, width=0.8, ax=ax13, alpha=alpha, yerr=[df.Spec_error_low_m, df.Spec_error_high_m], capsize=capsize)
    ax13.set_ylabel('Genes / Genome Size', fontweight='bold', fontsize=14)
    ax13.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax13.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax13.minorticks_on()
    ax13.set_axisbelow(True)
    ax13ytick = ax13.get_yticks()
    ax13.set_yticks(ax13ytick)

    df.Rare_norm.plot.bar(color=crare, width=0.8, ax=ax14, alpha=alpha, yerr=[df.Rare_error_low_m, df.Rare_error_high_m], capsize=capsize)
    ax14.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax14.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax14.minorticks_on()
    ax14.set_axisbelow(True)
    ax14.set_yticks(ax13ytick)

    df.Common_norm.plot.bar(color=ccommon, width=0.8, ax=ax15, alpha=alpha, yerr=[df.Common_error_low_m, df.Common_error_high_m], capsize=capsize)
    ax15.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax15.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax15.minorticks_on()
    ax15.set_axisbelow(True)
    ax15.set_yticks(ax13ytick)

    df.Core_norm.plot.bar(color=ccore, width=0.8, ax=ax16, alpha=alpha, yerr=[df.Core_error_low_m, df.Core_error_high_m], capsize=capsize)
    ax16.yaxis.grid(which="both", color=grid_color, linestyle='--', linewidth=1)
    #ax16.xaxis.grid(which="major", color=xgrid_color, linestyle='-', linewidth=1)
    ax16.minorticks_on()
    ax16.set_axisbelow(True)
    ax16.set_yticks(ax13ytick)

    ax13.set_xticklabels(df.index, rotation=30, fontsize=12, ha='right')
    ax14.set_xticklabels(df.index, rotation=30, fontsize=12, ha='right')
    ax15.set_xticklabels(df.index, rotation=30, fontsize=12, ha='right')
    ax16.set_xticklabels(df.index, rotation=30, fontsize=12, ha='right')

    plt.subplots_adjust(bottom=0.13, hspace=0.3, top=0.94, left=0.09, right=0.96)
    plt.savefig(f'{op}.png')
    plt.close()


def parse_pangenome_summary_file(cd, file, pangenome_summary):
    ''' parses a pangenome summary file, populates the dicts and returns them '''

    with open(f'{cd}{file}', 'r') as f:

        # Each pangenome summary file should only be one tab separated line
        X = f.readline().rstrip().split('\t')

        # Select the mean number of genes and pangenome for each organism
        p = float(X[1])
        m = float(X[6])
        # Make a list of the selected variables:
        # X[0] is Organism name. X[1] is Pangenome count, X[2] is core genome, X[3] common,
        # X[4] rare, X[5] specific, X[6] mean Genome Size (genes per genome), X[7] New genes (mean new genes per genome)
        # Also adds 95% Emperical Confidence Interval X[8] - X[19]
        pangenome_summary[X[0]] = [
            float(X[1]), float(X[1])/m, float(X[7]), float(X[7])/m,
            float(X[5])/p, float(X[4])/p, float(X[3])/p, float(X[2])/p,
            float(X[5])/m, float(X[4])/m, float(X[3])/m, float(X[2])/m,
            float(X[5]), float(X[4]), float(X[3]), float(X[2]),
            float(X[8]), float(X[9]), float(X[10]), float(X[11]),
            float(X[8])/m, float(X[9])/m, float(X[10])/m, float(X[11])/m,
            float(X[12]), float(X[13]), float(X[14]), float(X[15]),
            float(X[12])/p, float(X[13])/p, float(X[14])/p, float(X[15])/p,
            float(X[12])/m, float(X[13])/m, float(X[14])/m, float(X[15])/m,
            float(X[16]), float(X[17]), float(X[18]), float(X[19]),
            float(X[16])/p, float(X[17])/p, float(X[18])/p, float(X[19])/p,
            float(X[16])/m, float(X[17])/m, float(X[18])/m, float(X[19])/m
            ]
            
    return pangenome_summary


def compare_pangenomes(cd, op):
    ''' Reads the pangenome comparison files in cd, builds two dictionaries, and passes them to functions to build plots '''

    # Initialize Dictionary: order matches list l from parse_pangenome_summary_file() function.
    # frac is fraction of pangenome. Norm is normailized by mean genome size. Error is 95% empircal confidence interval.

    column_Order = [
        'Pangenome', 'Pan_norm', 'New_Genes', 'New_norm',
        'Specific_frac', 'Rare_frac', 'Common_frac', 'Core_frac',
        'Specific_norm', 'Rare_norm', 'Common_norm', 'Core_norm',
        'Specific', 'Rare', 'Common', 'Core',
        'Pan_error_low', 'Pan_error_high', 'New_error_low', 'New_error_high',
        'Pan_error_low_m', 'Pan_error_high_m', 'New_error_low_m', 'New_error_high_m',
        'Spec_error_low', 'Spec_error_high', 'Rare_error_low', 'Rare_error_high',
        'Spec_error_low_p', 'Spec_error_high_p', 'Rare_error_low_p', 'Rare_error_high_p',
        'Spec_error_low_m', 'Spec_error_high_m', 'Rare_error_low_m', 'Rare_error_high_m',
        'Common_error_low', 'Common_error_high', 'Core_error_low', 'Core_error_high',
        'Common_error_low_p', 'Common_error_high_p', 'Core_error_low_p', 'Core_error_high_p',
        'Common_error_low_m', 'Common_error_high_m', 'Core_error_low_m', 'Core_error_high_m'
        ]

    pangenome_summary = {'Order': column_Order}

    file_list = os.listdir(cd)

    # Read through files and populate dictionaries
    for file in file_list:
        if file == f'.DS_Store': continue
        pangenome_summary = parse_pangenome_summary_file(
            cd, file, pangenome_summary
            )

    # Pangenome Summary Plot
    df = pd.DataFrame(pangenome_summary).set_index(['Order']).T
    df = df.sort_values(by=['Pan_norm'], ascending=False)
    plot_pangenome_summary(df, op)
    print(df.index)

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
