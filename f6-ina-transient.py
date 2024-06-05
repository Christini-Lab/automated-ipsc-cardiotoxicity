import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from utility import get_valid_cells, get_vc_pts, get_model_response

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)


def plot_figure():
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6))
    fig.subplots_adjust(.1, .12, .95, .95)

    axs = axs.flatten()

    cell_name = 'flecainide_230403_008_1'
    plot_ex_cell(axs[0], channel='I_Na', cell_name=cell_name)
    ymin, ymax = axs[0].get_ylim()
    axs[0].annotate("Flec", xy=(2672.8, ymin + (ymax-ymin)*.5), color='k')
    axs[0].annotate("Wash", xy=(2677.7, ymin + (ymax-ymin)*.5), color='k')

    #cell_name = 'quinine_221206_003_2'
    #cell_name = 'quinine_221206_001_3'
    #cell_name = 'quinine_221206_001_3'
    cell_name = 'quinine_230403_011_1'
    plot_ex_cell(axs[1], channel='I_Na', cell_name=cell_name)
    ymin, ymax = axs[1].get_ylim()
    axs[1].annotate("Quin", xy=(2672.8, ymin + (ymax-ymin)*.5), color='k')
    axs[1].annotate("Wash", xy=(2677.7, ymin + (ymax-ymin)*.5), color='k')

    plot_curr(axs[2], 'flecainide', 'Na2_min', True)
    plot_curr(axs[3], 'quinine', 'Na2_min', True)

    letters = ['A', 'B', 'C', 'D']

    for i, ax in enumerate(axs):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(letters[i], y=.94, x=-.15)

    axs[2].set_ylabel(r"Change in $I_{Na}$ segment (%)")
    axs[2].set_ylim(-100, 120)
    axs[2].set_xlabel("Flecainide Concentration (xEFPC)")
    #axs[3].set_ylabel(r"Change in $I_{Na}$ segment (%)")
    axs[3].set_ylim(-100, 120)
    axs[3].set_xlabel("Quinine Concentration (xEFPC)")


    plt.savefig('figure-pdfs/f-ina-drug-block.pdf')

    plt.show()


def plot_ex_cell(ax, channel, cell_name):
    all_files = get_valid_cells(cell_name.split('_')[0]) 
    windows = {'I_Na': [2672, 2680],
               'I_Kr': [3670, 3730],
               'I_to': [4856, 4876],
               'I_CaL': [4245, 4260]}
    ch_ylims = {'I_Na': None,
                'I_Kr': [-20, 8],
                'I_to': [-30, 10], 
                'I_CaL': [-30, 5]}

    window = windows[channel]
    ylims = ch_ylims[channel]
    
    plot_window(cell_name, window, ax, ylims=ylims)

    ymin, ymax = ax.get_ylim()

    ax.annotate("", xy=(2674, ymax - (ymax-ymin)*.3), xytext=(2674, ymin+(ymax-ymin)*.1), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='red'))
    ax.annotate("", xy=(2677.5, ymin+(ymax-ymin)*.1), xytext=(2677.5, ymax - (ymax-ymin)*.3), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='grey'))
    #ax.annotate("", xy=(2674, -140), xytext=(2674, -10), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='grey'))

    ax.set_xlabel('Time (ms)')

    #plt.show()


def plot_curr(ax, drug_name, channel_name, is_change):
    all_file_names = {}
    all_file_names['dmso'] = get_valid_cells('dmso')
    all_file_names[drug_name] = get_valid_cells(drug_name)

    dmso_dat = []
    drug_dat = []

    for drug, all_f_names in all_file_names.items():
        for f_name in all_f_names:
            vc_summary = pd.read_csv(
                    f'./data/cells/{f_name}/vc_summary_stats.csv')

            if ('Na1' in channel_name) or ('Na2' in channel_name):
                if (vc_summary[vc_summary['compound'] == 'Baseline']
                                            [channel_name].mean() > -40):
                    continue

            new_list = [f_name, drug, channel_name]

            for comp_name in ['Baseline', '3xEFPC', '10xEFPC',
                                                '20xEFPC', 'Washout']:

                new_list.append(
                        vc_summary[vc_summary['compound'] == comp_name][channel_name].mean())

            if drug == drug_name:
                drug_dat.append(new_list)
            else:
                dmso_dat.append(new_list)

    dmso_dat = pd.DataFrame(dmso_dat, columns=['file', 'compound', 'feature', 'baseline', '3x', '10x', '20x', 'wash'])
    drug_dat =  pd.DataFrame(drug_dat, columns=['file', 'compound', 'feature', 'baseline', '3x', '10x', '20x', 'wash'])

    if is_change:
        new_drug_cols = [100*(drug_dat.iloc[:, 3+i] - drug_dat['baseline'])/drug_dat['baseline'] for i in range(0, 5)]
        drug_dat['baseline'], drug_dat['3x'], drug_dat['10x'], drug_dat['20x'], drug_dat['wash'] = new_drug_cols[0], new_drug_cols[1], new_drug_cols[2], new_drug_cols[3], new_drug_cols[4]

        new_dmso_cols = [100*(dmso_dat.iloc[:, 3+i] - dmso_dat['baseline'])/dmso_dat['baseline'] for i in range(0, 5)]
        dmso_dat['baseline'], dmso_dat['3x'], dmso_dat['10x'], dmso_dat['20x'], dmso_dat['wash'] = new_dmso_cols[0], new_dmso_cols[1], new_dmso_cols[2], new_dmso_cols[3], new_dmso_cols[4]
    else:
        new_drug_cols = [drug_dat.iloc[:, 3+i] - drug_dat['baseline'] for i in range(0, 5)]
        drug_dat['baseline'], drug_dat['3x'], drug_dat['10x'], drug_dat['20x'], drug_dat['wash'] = new_drug_cols[0], new_drug_cols[1], new_drug_cols[2], new_drug_cols[3], new_drug_cols[4]

        new_dmso_cols = [dmso_dat.iloc[:, 3+i] - dmso_dat['baseline'] for i in range(0, 5)]
        dmso_dat['baseline'], dmso_dat['3x'], dmso_dat['10x'], dmso_dat['20x'], dmso_dat['wash'] = new_dmso_cols[0], new_dmso_cols[1], new_dmso_cols[2], new_dmso_cols[3], new_dmso_cols[4]

    x_arr = np.array([0, 1, 2, 3, 4])



    #FIX THIS DMSO_DIFF ISSUE
    ax.errorbar(x_arr, dmso_dat.iloc[:, 3:].mean(0), dmso_dat.iloc[:, 3:].std(0), linestyle='-', c='k', label='DMSO', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='k', alpha=.4) for y_vals in dmso_dat.iloc[:, 3:].values]

    x_arr = x_arr + 0.05
    ax.errorbar(x_arr, drug_dat.iloc[:, 3:].mean(0), drug_dat.iloc[:, 3:].std(0), linestyle='--', c='r', label=drug.capitalize(), capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='r', alpha=.4) for y_vals in drug_dat.iloc[:, 3:].values]

    max_val = np.concatenate((dmso_dat.iloc[:, 3:].values,
                                drug_dat.iloc[:, 3:].values)).max()
    min_val = np.concatenate((dmso_dat.iloc[:, 3:].values,
                                    drug_dat.iloc[:, 3:].values)).min()

    y_sig = max_val + np.abs(max_val)*.1
    is_sig = [ttest_ind(dmso_dat.iloc[:, 3:].values[:, i],
                        drug_dat.iloc[:, 3:].values[:, i]).pvalue < .05
                            for i in range(0, 5)]

    #[ax.text(i, y_sig, '*', fontsize=14) for i in range(0, 5) if is_sig[i]]

    sig_vals = [ttest_ind(dmso_dat.iloc[:, 3:].values[:, i],
                    drug_dat.iloc[:, 3:].values[:, i]).pvalue
                        for i in range(0, 4)]


    for i, sig_val in enumerate(sig_vals):
        if sig_val < .001:
            aster = '***'
        elif sig_val < .01:
            aster = '**'
        elif sig_val < .05:
            aster = '*'
        else:
            aster = ''

        ax.text(i-.05, y_sig, aster, fontsize=14)

    #ax.set_ylabel(r'$I_{Na}$ Change from Baseline (%)')

    #ax.legend()

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(['B', '3x', '10x', '20x', 'W'])
    #ax.set_ylim(min_val - np.abs(min_val)*.1, max_val + np.abs(max_val*.25))


# PLOTTING FUNCTIONS
def plot_window(f, window, ax_c, ylims=None):
    
    idx = np.array([int(val*25) for val in window])

    vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
    vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')

    vc_slice = vc_dat.iloc[idx[0]:idx[1], :]

    times = np.linspace(0, int(vc_dat.shape[0]/25), vc_dat.shape[0]+1)[0:-1]
    times = times[idx[0]:idx[1]]

    for i, k in enumerate(vc_dat.keys()):
        curr_dat = vc_slice[k].values

        curr_meta = vc_meta[vc_meta['sweep'] == k]

        if 'Baseline' == curr_meta['compound'].values[0]:
            cols = c=(0, 0, 0)
        if '3xEFPC' == curr_meta['compound'].values[0]:
            cols = c=(1, .7, .7)
        if '10xEFPC' == curr_meta['compound'].values[0]:
            cols =(1, .4, .4)
        if '20xEFPC' == curr_meta['compound'].values[0]:
            cols =(1, 0, 0)
        if 'Washout' == curr_meta['compound'].values[0]:
            cols =(.5, .5, .5)

        if i == 0:
            ax_c.plot(times, curr_dat/vc_meta['cm'].mean(), c=cols)
        else:
            ax_c.plot(times, curr_dat/vc_meta['cm'].mean(), c=cols)

    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)
    ax_c.set_ylabel('A/F')
    try:
        ax_c.set_xlim(times[0], times[-1])
    except:
        import pdb
        pdb.set_trace()
    if ylims is not None:
        ax_c.set_ylim(ylims[0], ylims[1])

    scales = {'Baseline': {'g_Kr': 1, 'g_to': 1, 'g_Na': 1, 'g_CaL': 1},
              '3xEFPC': {'g_Kr': .3, 'g_to': .75, 'g_Na': .9, 'g_CaL': .98},
              '10xEFPC': {'g_Kr': .85, 'g_to': .45, 'g_Na': .52, 'g_CaL': .87},
              '20xEFPC': {'g_Kr': .9, 'g_to': .6, 'g_Na': .8, 'g_CaL': .3}
              }

    #for conc, scale_vals in scales.items():
    #    dat = get_model_response('Kernik', scale_vals)


def main():
    plot_figure()


if __name__ == '__main__':
    main()
