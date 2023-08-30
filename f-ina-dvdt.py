import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import pointplot, swarmplot 
from scipy.stats import ttest_ind


from utility import get_valid_cells
import matplotlib

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)
matplotlib.rcParams['pdf.fonttype'] = 42



def plot_figure_dvdt():
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6))
    fig.subplots_adjust(.1, .15, .95, .95)

    axs = axs.flatten()


    plot_dvdt(axs[0], drug_name='flecainide', is_pct_change=True)
    plot_dvdt(axs[1], drug_name='quinine', is_pct_change=True)

    window = [2674.3, 2685]
    channel = 'I_Na'

    plot_curr(axs[2], 'flecainide', window, channel, is_change=True)
    plot_curr(axs[3], 'quinine', window, channel, is_change=True)

    for vals in [[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)
    
    axs[0].text(1.15, 270, 'Flecainide', fontsize=14)
    axs[1].text(1.15, 270, 'Quinine', fontsize=14)

    axs[0].set_ylabel(r'$dV/dt_{max}$')

    plt.savefig('./figure-pdfs/f-flec_quin-dvdt-ina.pdf')
    plt.show()


def plot_dvdt(ax, drug_name, is_pct_change=False):
    dmso_means = []
    drug_means = []

    dmso_meta = get_all_ap_meta('dmso')
    drug_meta = get_all_ap_meta(drug_name)

    for meta in [dmso_meta, drug_meta]:
        for f, ap_meta in meta.items():
            if ap_meta['is_upstroke'].sum() < 45:
                continue

            mean_dvdt = []
            for compound in ap_meta['compound'].unique():
                mean_dvdt.append([ap_meta[ap_meta[
                                    'compound'] == compound]['dvdt_max'].mean(),
                                  ap_meta[ap_meta[
                                    'compound'] == compound]['dvdt_max'].std()])

            mean_dvdt = np.array(mean_dvdt)
            if is_pct_change:
                curr_dvdt = (mean_dvdt[:, 0] - mean_dvdt[0,0])/mean_dvdt[0,0] * 100
            else:
                curr_dvdt = mean_dvdt[:, 0]


            if 'dmso' in f:
                dmso_means.append(curr_dvdt)
            else:
                drug_means.append(curr_dvdt)

            print(f)

    x_arr = np.array([0, 1, 2, 3, 4])
    ax.errorbar(x_arr-.05, np.array(dmso_means).mean(0), np.array(dmso_means).std(0), linestyle='-', c='k', label='DMSO', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08)-.05, y_vals, c='k', alpha=.4) for y_vals in dmso_means]
    ax.errorbar(x_arr+0.05, np.array(drug_means).mean(0), np.array(drug_means).std(0), linestyle='--', c='r', label=drug_name.capitalize(), capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08)+.05, y_vals, c='r', alpha=.4) for y_vals in drug_means]

    #ax.set_ylabel(r'$dV/dt_{max}$ (V/s)')
    dmso_means = np.array(dmso_means)
    drug_means = np.array(drug_means)

    max_val = np.concatenate((dmso_means, drug_means)).max()
    min_val = np.concatenate((dmso_means, drug_means)).min()

    y_sig = max_val + np.abs(max_val)*.1
    is_sig = [ttest_ind(dmso_means[:,i], drug_means[:,i]).pvalue < .05 for i in range(0, 4)]
    sig_vals = [ttest_ind(dmso_means[:,i], drug_means[:,i]).pvalue for i in range(0, 4)]

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
    
    #[ax.text(i, y_sig, '*', fontsize=14) for i in range(0, 4) if is_sig[i]]


    #ax.legend()
    ax.set_title(drug_name.capitalize())
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Wash'])

    ax.set_ylim(min_val - np.abs(min_val)*.1, max_val + np.abs(max_val*.25))


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

        if channel == 'I_Na':
            if 'dmso' in f:
                temp_window = [1550.5, 1580]
                idx = np.array([int(val*25) for val in temp_window])
            else:
                idx = np.array([int(val*25) for val in window])

        vc_slice = vc_dat.iloc[idx[0]:idx[1], :]

        min_vals = vc_slice.min().values / vc_meta['cm'].values
        min_vals = [np.mean(min_vals[i*2:i*2+2]) for i in range(0, 5)]

        if channel == 'I_Na':
            if np.min(min_vals) > -20:
                print(f'{f}: No sign of Na at any concentration')
                continue

        if 'dmso' in f:
            vc_dict['Drug'].append('DMSO')
        else:
            vc_dict['Drug'].append(drug)

        for i, compound in enumerate(['Baseline', '3xEFPC', '10xEFPC',
                                                  '20xEFPC', 'Washout']):
            vc_dict[compound].append(min_vals[i])

        vc_dict['File'].append(f)

        print(f)

    dat = pd.DataFrame(vc_dict)

    dmso_dat = dat[dat['Drug'] == 'DMSO']
    drug_dat = dat[dat['Drug'] == drug] 

    dmso_diff = 100*np.array([((vals-vals[0])/vals[0]).astype('float64') for vals in dmso_dat.values[:, 2:]])
    #dmso_diff = np.diff(dmso_dat.values[:, 2:]).astype('float64')
    x_arr = np.array([0, 1, 2, 3, 4])

    #FIX THIS DMSO_DIFF ISSUE
    ax.errorbar(x_arr, dmso_diff[:-1,:].mean(0), dmso_diff[:-1,:].std(0), linestyle='-', c='k', label='DMSO', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='k', alpha=.4) for y_vals in dmso_diff]

    drug_diff = np.diff(drug_dat.values[:, 2:]).astype('float64')
    drug_diff = 100*np.array([((vals-vals[0])/vals[0]).astype('float64') for vals in drug_dat.values[:, 2:]])
    x_arr = x_arr + 0.05
    ax.errorbar(x_arr, drug_diff.mean(0), drug_diff.std(0), linestyle='--', c='r', label=drug.capitalize(), capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='r', alpha=.4) for y_vals in drug_diff]

    max_val = np.concatenate((dmso_diff, drug_diff)).max()
    min_val = np.concatenate((dmso_diff, drug_diff)).min()

    y_sig = max_val + np.abs(max_val)*.1
    is_sig = [ttest_ind(dmso_diff[:,i], drug_diff[:,i]).pvalue < .05 for i in range(0, 5)]

    #[ax.text(i, y_sig, '*', fontsize=14) for i in range(0, 5) if is_sig[i]]

    sig_vals = [ttest_ind(dmso_diff[:,i], drug_diff[:,i]).pvalue for i in range(0, 4)]

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

    ax.set_ylabel(r'Change from Baseline (%)')

    #ax.legend()

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(['B', '3x', '10x', '20x', 'W'])
    ax.set_ylim(min_val - np.abs(min_val)*.1, max_val + np.abs(max_val*.25))


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


def get_all_ap_meta(drug_name):
    all_files = get_valid_cells(drug_name)

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

