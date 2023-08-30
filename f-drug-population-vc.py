import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utility import get_valid_cells, get_vc_pts

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)


def plot_drug_channel_change(drug, channel):
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3))
    fig.subplots_adjust(.1, .1, .95, .9, wspace=.3)
    axs = axs.flatten()
    windows = {'I_Na': [1556.4, 1580],
               'I_Kr': [3670, 3680],
               'I_to': [4862, 4872]}

    window = windows[channel]

    plot_curr(axs[0], drug, window, channel)
    plot_curr(axs[1], drug, window, channel, is_change=True)

    for vals in [[0, 'A'], [1, 'B']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    fig.suptitle(f'{drug} effect on {channel}')

    plt.savefig(f'./figure-pdfs/f-dmso-vs-{drug}-vc{channel}.pdf')
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

        vc_slice = vc_dat.iloc[idx[0]:idx[1], :]

        min_vals = vc_slice.min().values / vc_meta['cm'].values
        min_vals = [np.mean(min_vals[i*2:i*2+2]) for i in range(0, 5)]

        if channel == 'I_Na':
            if np.min(min_vals) > -20:
                print(f'{f}: No sign of Na at any concentration')
                continue

        if 'dmso' in f:
            vc_dict['Drug'].append('DMSO')
        if 'flecainide' in f:
            vc_dict['Drug'].append('Flecainide')
    
        for i, compound in enumerate(['Baseline', '3xEFPC', '10xEFPC',
                                                  '20xEFPC', 'Washout']):
            vc_dict[compound].append(min_vals[i])

        vc_dict['File'].append(f)

        print(f)

    dat = pd.DataFrame(vc_dict)

    dmso_dat = dat[dat['Drug'] == 'DMSO']
    flec_dat = dat[dat['Drug'] == 'Flecainide']

    if is_change:
        dmso_diff = 100*np.array([((vals-vals[0])/vals[0]).astype('float64') for vals in dmso_dat.values[:, 2:]])
        #dmso_diff = np.diff(dmso_dat.values[:, 2:]).astype('float64')
        x_arr = np.array([0, 1, 2, 3, 4])
        #FIX THIS DMSO_DIFF ISSUE
        ax.errorbar(x_arr, dmso_diff[:-1,:].mean(0), dmso_diff[:-1,:].std(0), linestyle='-', c='k', label='DMSO', capsize=3)
        [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='k', alpha=.4) for y_vals in dmso_diff]

        flec_diff = np.diff(flec_dat.values[:, 2:]).astype('float64')
        flec_diff = 100*np.array([((vals-vals[0])/vals[0]).astype('float64') for vals in flec_dat.values[:, 2:]])
        x_arr = x_arr + 0.05
        ax.errorbar(x_arr, flec_diff.mean(0), flec_diff.std(0), linestyle='-', c='r', label='Flecainide', capsize=3)
        [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='r', alpha=.4) for y_vals in flec_diff]

        ax.set_ylabel(r'Change from Baseline (%)')

        ax.legend()

        #ax.set_xticks([0, 1, 2, 3, 4])
        #ax.set_xticklabels(['B->3x', '3x->10x', '10x-20x', '20x->W'])
        ax.set_xticks([0, 1, 2, 3, 4])
        ax.set_xticklabels(['B', '3x', '10x', '20x', 'W'])
    else:
        x_arr = np.array([0, 1, 2, 3, 4])
        ax.errorbar(x_arr, dmso_dat.mean().values, dmso_dat.std().values, linestyle='-', c='k', label='DMSO', capsize=3)
        [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='k', alpha=.4) for y_vals in dmso_dat.values[:, 2:]]
        x_arr = x_arr + 0.05
        ax.errorbar(x_arr, flec_dat.mean().values, flec_dat.std().values, linestyle='-', c='r', label='Flecainide', capsize=3)
        [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='r', alpha=.4) for y_vals in flec_dat.values[:, 2:]]

        ax.set_ylabel(r'Current (A/F)')

        ax.legend()

        ax.set_xticks([0, 1, 2, 3, 4])
        ax.set_xticklabels(['B', '3x', '10x', '20x', 'W'])



# PLOTTING FUNCTIONS
def plot_window(f, window, ax_c, ax_v=None):
    if ax_v is not None:
        vc_pts = get_vc_pts()
        ax_v.plot(vc_pts[:, 0], vc_pts[:, 1], 'k')
        ax_v.set_ylabel('mV')
        ax_v.spines['top'].set_visible(False)
        ax_v.spines['right'].set_visible(False)
        ax_v.set_ylim(-50, 50)
    
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
            ax_c.plot(times, curr_dat/vc_meta['cm'].mean(), c=cols, label=f)
        else:
            ax_c.plot(times, curr_dat/vc_meta['cm'].mean(), c=cols)

    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)
    ax_c.set_ylabel('A/F')
    ax_c.set_xlim(times[0], times[-1])
    ax_c.set_ylim(-20, 8)
    ax_c.legend(loc=4)



def main():
    plot_drug_channel_change('flecainide', 'I_Na')
    #plot_drug_channel_change('flecainide', 'I_Kr')
    #plot_drug_channel_change('flecainide', 'I_to')
    


if __name__ == '__main__':
    main()

