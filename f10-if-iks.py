import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import pointplot, swarmplot 
import pickle
from scipy.stats import ttest_ind

from utility import get_valid_cells

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)


def plot_figure():
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6))
    fig.subplots_adjust(.1, .15, .95, .95, wspace=.3)

    axs = axs.flatten()

    channel = 'I_f'
    window = [6850, 6885]
    is_change = True #Abs
    plot_curr(axs[0], 'flecainide', window, channel, is_change=is_change)
    plot_curr(axs[1], 'quinine', window, channel, is_change=is_change)


    channel = 'I_Ks'
    window = [9350, 9400]
    is_change = True #Abs
    plot_curr(axs[2], 'flecainide', window, channel, is_change=is_change)
    plot_curr(axs[3], 'quinine', window, channel, is_change=is_change)

    for vals in [[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    axs[0].text(1.15, 8, 'Flecainide', fontsize=14)
    axs[0].set_ylim(-25, 10)

    axs[1].text(1.15, 8, 'Quinine', fontsize=14)
    axs[1].set_ylim(-20, 10)


    axs[0].set_ylabel(r'Change in $I_{f}$ segment (A/F)')
    axs[2].set_ylabel(r'Change in $I_{Ks}$ segment (A/F)')

    plt.savefig(f'./figure-pdfs/f-flec_quin-I_f-I_Ks.pdf')
    plt.show()


def plot_curr(ax, drug, window, channel, is_change=False):
    all_files = get_valid_cells('dmso') + get_valid_cells(drug)

    idx = np.array([int(val*25) for val in window])

    vc_dict = {'File': [],
               'Drug': [],
               'Baseline': [],
               '3xEFPC': [],
               '10xEFPC': [],
               '20xEFPC': [],
               'Washout': []}

    for f in all_files:
        vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
        vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')

        idx = np.array([int(val*25) for val in window])

        vc_slice = vc_dat.iloc[idx[0]:idx[1], :]

        mean_vals = vc_slice.mean().values / vc_meta['cm'].values
        mean_vals = [np.mean(mean_vals[i*2:i*2+2]) for i in range(0, 5)]

        if 'dmso' in f:
            vc_dict['Drug'].append('DMSO')
        else:
            vc_dict['Drug'].append(drug)

        for i, compound in enumerate(['Baseline', '3xEFPC', '10xEFPC',
                                                  '20xEFPC', 'Washout']):
            vc_dict[compound].append(mean_vals[i])

        vc_dict['File'].append(f)

        print(f)

    dat = pd.DataFrame(vc_dict)

    dmso_dat = dat[dat['Drug'] == 'DMSO']
    drug_dat = dat[dat['Drug'] == drug] 

    if is_change:
        dmso_diff = np.array([(vals-vals[0]).astype('float64') for vals in dmso_dat.values[:, 2:]])
        drug_diff = np.array([(vals-vals[0]).astype('float64') for vals in drug_dat.values[:, 2:] if vals[2]-vals[0] <15])
        #ax.set_ylabel(r'Change from Baseline (A/F)')
    else:
        dmso_diff = np.array([100*(vals-vals[0]).astype('float64')/vals[0] for vals in dmso_dat.values[:, 2:]])
        drug_diff = np.array([100*(vals-vals[0]).astype('float64')/vals[0] for vals in drug_dat.values[:, 2:] if vals[2]-vals[0] <15])
        #ax.set_ylabel(r'Change from Baseline (%)')

    #dmso_diff = np.diff(dmso_dat.values[:, 2:]).astype('float64')
    x_arr = np.array([0, 1, 2, 3, 4])

    ax.errorbar(x_arr, dmso_diff[:-1,:].mean(0), dmso_diff[:-1,:].std(0), linestyle='-', c='k', label='DMSO', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='k', alpha=.4) for y_vals in dmso_diff]

    #drug_diff = np.diff(drug_dat.values[:, 2:]).astype('float64')
    x_arr = x_arr + 0.05
    ax.errorbar(x_arr, drug_diff.mean(0), drug_diff.std(0), linestyle='--', c='r', label=drug.capitalize(), capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='r', alpha=.4) for y_vals in drug_diff]


    max_val = np.concatenate((dmso_diff, drug_diff)).max()
    min_val = np.concatenate((dmso_diff, drug_diff)).min()

    y_sig = max_val + np.abs(max_val)*.1
    is_sig = [ttest_ind(dmso_diff[:,i], drug_diff[:,i]).pvalue < .05 for i in range(0, 5)]

    [ax.text(i, y_sig, '*', fontsize=18) for i in range(0, 5) if is_sig[i]]


    #ax.legend()

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(['B', '3x', '10x', '20x', 'W'])
    ax.set_ylim(min_val - np.abs(min_val)*.1, max_val + np.abs(max_val*.25))


# Utility
def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]
    return np.array(new_vals)


def main():
    plot_figure()


if __name__ == '__main__':
    main()

