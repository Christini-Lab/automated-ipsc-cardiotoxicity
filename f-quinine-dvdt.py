import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import pointplot, swarmplot 

from utility import get_valid_cells


plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)


def plot_figure_dvdt():
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3))
    fig.subplots_adjust(.1, .15, .95, .95)

    axs = axs.flatten()

    meta = get_all_ap_meta()

    plot_dvdt(axs[0], meta)

    plot_dvdt_changes(axs[1], meta)

    for vals in [[0, 'A'], [1, 'B']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    plt.savefig('./figure-pdfs/f-dvdt-quinine.pdf')
    plt.show()


def plot_dvdt(ax, meta):
    dmso_means = []
    flecainide_means = []

    for f, ap_meta in meta.items():
        if ap_meta['is_upstroke'].sum() < 45:
            continue

        mean_dvdt = []
        for compound in ap_meta['compound'].unique():
            mean_dvdt.append([ap_meta[ap_meta['compound'] == compound]['dvdt_max'].mean(),
                               ap_meta[ap_meta['compound'] == compound]['dvdt_max'].std()])

        mean_dvdt = np.array(mean_dvdt)

        if 'dmso' in f:
            dmso_means.append(mean_dvdt[:, 0])
        else:
            flecainide_means.append(mean_dvdt[:, 0])

        print(f)

    x_arr = np.array([0, 1, 2, 3, 4])
    ax.errorbar(x_arr-.05, np.array(dmso_means).mean(0), np.array(dmso_means).std(0), linestyle='-', c='k', label='DMSO', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08)-.05, y_vals, c='k', alpha=.4) for y_vals in dmso_means]
    ax.errorbar(x_arr+0.05, np.array(flecainide_means).mean(0), np.array(flecainide_means).std(0), linestyle='--', c='r', label='Flecainide', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08)+.05, y_vals, c='r', alpha=.4) for y_vals in flecainide_means]


    ax.set_ylabel(r'$dV/dt_{max}$ (V/s)')

    ax.legend()

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Wash'])


def plot_dvdt_changes(ax, meta):
    """
    TODO:
        Currently, I exclude APs that do not generate an upstroke.
        I SHOULD exclude APs if there is no noticable INa in the VC data
    """
    mean_dmso_changes = []
    mean_flec_changes = []
    
    trial_changes_flec = {'Baseline->3xEFPC': [],
                     '3xEFPC->10xEFPC': [],
                     '10xEFPC->20xEFPC': [],
                     '20xEFPC->Washout': []} 
    trial_changes_dmso = {'Baseline->3xEFPC': [],
                     '3xEFPC->10xEFPC': [],
                     '10xEFPC->20xEFPC': [],
                     '20xEFPC->Washout': []} 
    trial_changes = {'Drug': [],
                     'Concentration': [],
                     'dvdt_change': []}

    for f, ap_meta in meta.items():
        if ap_meta['is_upstroke'].sum() < 45:
            continue

        split_sweeps = [x.split('_') for x in ap_meta['sweep']]
        arr = np.sort(np.unique(np.array(split_sweeps)[:, 1]).astype(int))

        for i, trial_num in enumerate(arr):
            curr_set = ap_meta[ap_meta['sweep'].str.contains(f'1_{trial_num}')]
            compounds = curr_set['compound'].unique()

            mean_start = curr_set[curr_set['compound'] == compounds[0]]['dvdt_max'].mean()
            mean_end = curr_set[curr_set['compound'] == compounds[1]]['dvdt_max'].mean()

            std_start = curr_set[curr_set['compound'] == compounds[0]]['dvdt_max'].std()
            std_end = curr_set[curr_set['compound'] == compounds[1]]['dvdt_max'].std()

            diff = mean_end - mean_start
            diff_std = np.sqrt(std_start**2+std_end**2)

            #if compounds[1] == 'Washout':
            #    print(f'{f}: {diff}')
            #    if diff < -20:
            #        import pdb
            #        pdb.set_trace()

            trial_changes['Concentration'].append(f'{compounds[0]}->{compounds[1]}')
            trial_changes['dvdt_change'].append(diff)

            if 'dmso' in f:
                trial_changes_dmso[f'{compounds[0]}->{compounds[1]}'].append(diff)
                trial_changes['Drug'].append('DMSO')
            else:
                trial_changes_flec[f'{compounds[0]}->{compounds[1]}'].append(diff)
                trial_changes['Drug'].append('Flecainide')

    dmso_df = pd.DataFrame(trial_changes_dmso).melt()
    flec_df = pd.DataFrame(trial_changes_flec).melt()
    total_changes = pd.DataFrame(trial_changes)

    g = swarmplot(x='Concentration', y='dvdt_change', hue='Drug', data=total_changes, palette=['k', 'r'], ax=ax, zorder=1, alpha=.4)
    g.legend_.remove()
    g = pointplot(x='Concentration', y='dvdt_change', hue='Drug', data=total_changes, join=False, capsize=.05, markers='_', ax=ax, palette=['k', 'r'], ci='sd', dodge=True)
    #g._legend.remove()

    #swarmplot(x='variable', y='value', data=dmso_df, color='k', ax=ax, zorder=1, alpha=.4)
    #pointplot(x='variable', y='value', data=dmso_df, join=False, capsize=.05, markers='_', ax=ax, color='k', ci='sd')

    #swarmplot(x='variable', y='value', data=flec_df, color='r', ax=ax, zorder=1, alpha=.4)
    #pointplot(x='variable', y='value', data=flec_df, join=False, capsize=.05, markers='_', ax=ax, color='r', ci='sd')

    ax.set_ylabel(r'$\Delta$$dV/dt_{max}$ (V/s)')

    #ax.legend()
    #for i in [0, 1]:
    #    axs[i].set_xticks([0, 1, 2, 3, 4])
    #    axs[i].set_xticklabels(['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Wash'])
    #for i in [2, 3]:
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(['B->3x', '3x->10x', '10x->20x', '20x->W'])
    ax.legend().legendHandles[0].set_label('')
    ax.legend().legendHandles[1].set_label('')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[2:], ['DMSO', 'Flecainide'])



# Utility
def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]
    return np.array(new_vals)


def get_dvdt(ap_num, ap):
    stim_region = ap.values[int(48*25):int(65*25)]

    vals = moving_average(stim_region, 7)
    times = np.linspace(48,65, len(stim_region))
    times_avg = np.linspace(48, 65, len(vals))

    diff = np.diff(vals) / (times_avg[2]-times_avg[1])

    max_dvdt = np.max(diff)

    stim_start_idx = np.argmin(np.abs(times_avg - 50.5))
    stim_start_dvdt = np.mean(diff[stim_start_idx:(stim_start_idx+3)])

    max_dvdt -= stim_start_dvdt

    is_na_during_stim = not np.all(np.abs(np.diff(diff[7:20])) < 10)
    is_late_upstroke = not np.max(diff[21:35]) < 5

    if is_na_during_stim or is_late_upstroke:
        is_na = True
    else:
        is_na = False

    #fig, axs = plt.subplots(3, 1, sharex=True)
    #axs[0].plot(times, stim_region)
    #axs[0].scatter(times_avg, vals)
    #axs[1].plot(times_avg[1:], diff/(times_avg[2]-times_avg[1]))
    #axs[2].plot(times_avg[2:], np.diff(diff))
    #plt.show()

    #import pdb
    #pdb.set_trace()

    return [max_dvdt, is_na]


def get_all_ap_meta():
    dmso_files = get_valid_cells('dmso')
    flecainide_files = get_valid_cells('quinine')

    all_files = dmso_files + flecainide_files
    #all_files = flecainide_files

    meta = {}

    
    for f in all_files:
        #f = 'dmso_220829_001_2' #NO NaV
        ap_dat = pd.read_csv(f'./data/cells/{f}/ap_df.csv')
        ap_meta = pd.read_csv(f'./data/cells/{f}/ap_meta.csv')

        dvdt_is_na = np.array([get_dvdt(k, ap_dat[k]) for k in ap_meta['sweep'].values])
        ap_meta['dvdt_max'] = dvdt_is_na[:, 0]
        ap_meta['is_upstroke'] = dvdt_is_na[:, 1]

        meta[f] = ap_meta

        print(f'\t There are {ap_meta["is_upstroke"].sum()} upstrokes.')
        print(dvdt_is_na[:, 0])

    return meta


def main():
    plot_figure_dvdt()


if __name__ == '__main__':
    main()
