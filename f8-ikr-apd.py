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


def plot_figure_apd():
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6))
    fig.subplots_adjust(.1, .15, .95, .95, wspace=.3)

    axs = axs.flatten()

    plot_apd(axs[0], drug_name='flecainide')
    plot_apd(axs[1], drug_name='quinine')

    channel = 'I_Kr'
    window = [3690, 3695]
    is_change = True #Abs

    #channel = 'I_Ks'
    #window = [9200, 9300]
    #is_change = False #Percent

    #channel = 'I_to'
    #window = [4862, 4880]
    #is_change = False #Percent


    plot_curr(axs[2], 'flecainide', window, channel, is_change=is_change)
    plot_curr(axs[3], 'quinine', window, channel, is_change=is_change)

    for vals in [[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    axs[0].text(1.15, 62, 'Flecainide', fontsize=14)
    axs[1].text(1.15, 120, 'Quinine', fontsize=14)


    axs[0].set_ylabel(r'Change in $APD_{90}$ (%)')

    plt.savefig(f'./figure-pdfs/f-flec_quin-apd-{channel}.pdf')
    plt.show()


def plot_apd(ax, drug_name):
    dmso_means = []
    drug_means = []

    dmso_meta = pickle.load(open('./data/dmso_apd_meta.pkl', 'rb'))
    #get_all_ap_meta('dmso')
    drug_meta = pickle.load(open(f'./data/{drug_name}_apd_meta.pkl', 'rb'))
    #get_all_ap_meta(drug_name)

    for meta in [dmso_meta, drug_meta]:
        for f, ap_meta in meta.items():
            if ap_meta['is_upstroke'].sum() < 45:
                continue

            mean_dvdt = []
            for compound in ap_meta['compound'].unique():
                mean_dvdt.append([ap_meta[ap_meta[
                                    'compound'] == compound]['apd90'].mean(),
                                  ap_meta[ap_meta[
                                    'compound'] == compound]['apd90'].std()])

            mean_dvdt = np.array(mean_dvdt)
            to_append = 100*(mean_dvdt[:, 0] - mean_dvdt[0][0]) / mean_dvdt[0][0]

            if 'dmso' in f:
                #dmso_means.append(mean_dvdt[:, 0])
                dmso_means.append(to_append)
            else:
                #drug_means.append(mean_dvdt[:, 0])
                drug_means.append(to_append)

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

    [ax.text(i, y_sig, '*', fontsize=18) for i in range(0, 4) if is_sig[i]]

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

        idx = np.array([int(val*25) for val in window])

        vc_slice = vc_dat.iloc[idx[0]:idx[1], :]

        min_vals = vc_slice.min().values / vc_meta['cm'].values
        min_vals = [np.mean(min_vals[i*2:i*2+2]) for i in range(0, 5)]

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

    if is_change:
        dmso_diff = np.array([(vals-vals[0]).astype('float64') for vals in dmso_dat.values[:, 2:]])
        drug_diff = np.array([(vals-vals[0]).astype('float64') for vals in drug_dat.values[:, 2:] if vals[2]-vals[0] <15])
        ax.set_ylabel(r'Change in $I_{Kr}$ segment (A/F)')
    else:
        dmso_diff = np.array([100*(vals-vals[0]).astype('float64')/vals[0] for vals in dmso_dat.values[:, 2:]])
        drug_diff = np.array([100*(vals-vals[0]).astype('float64')/vals[0] for vals in drug_dat.values[:, 2:] if vals[2]-vals[0] <15])
        ax.set_ylabel(r'Change from Baseline (%)')

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
    sig_vals = [ttest_ind(dmso_diff[:,i], drug_diff[:,i]).pvalue for i in range(0, 5)]

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


def get_apd(ap_num, ap):
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

    ap_orig = ap[:11000]
    times_orig = np.linspace(0, len(ap_orig)/25, len(ap_orig))

    ap = moving_average(ap[:11000], 4)

    max_v = np.max(ap) 
    apa = max_v - np.min(ap)
    apd90_v = max_v - apa*.9
    argmax_v = np.argmax(ap)

    apd90_idx = np.argmin(np.abs(ap[argmax_v:] - apd90_v)) + argmax_v
    dvdt_max = np.argmax(np.diff(ap))

    #plt.plot(ap)
    #plt.axvline(dvdt_max, color='r')
    #plt.axvline(apd90_idx, color='b')

    times_new = np.linspace(0, 11000/25, len(ap))
    apd90 = times_new[apd90_idx] - times_new[dvdt_max]

    return [apd90, is_na]


def get_all_ap_meta(drug_name):
    all_files = get_valid_cells(drug_name)

    meta = {}
    
    for f in all_files:
        #f = 'dmso_220829_001_2' #NO NaV
        ap_dat = pd.read_csv(f'./data/cells/{f}/ap_df.csv')
        ap_meta = pd.read_csv(f'./data/cells/{f}/ap_meta.csv')

        apd_is_na = np.array([get_apd(k, ap_dat[k]) for k in ap_meta['sweep'].values])
        ap_meta['apd90'] = apd_is_na[:, 0]
        ap_meta['is_upstroke'] = apd_is_na[:, 1]

        meta[f] = ap_meta

        print(f'\t There are {ap_meta["is_upstroke"].sum()} upstrokes.')
        print(apd_is_na[:, 0])

    pickle.dump(meta, open(f'./data/{drug_name}_apd_meta.pkl', 'wb'))

    return meta



def main():
    plot_figure_apd()


if __name__ == '__main__':
    main()

