import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import regplot
from scipy import stats
import matplotlib

from utility import get_valid_cells


plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)
matplotlib.rcParams['pdf.fonttype'] = 42



def plot_single_dvdt_vs_vc():
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.5))
    fig.subplots_adjust(.1, .15, .95, .95)

    plot_dvdt_vs_vc(axs[0], drug_name='flecainide')
    plot_dvdt_vs_vc(axs[1], drug_name='quinine')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    for vals in [[0, 'A'], [1, 'B']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)


    axs[0].set_ylabel(r'$dV/dt_{max}$ (V/s)')
    axs[0].set_xlabel(r'Current (A/F)')
    axs[1].set_xlabel(r'Current (A/F)')
    #ax.set_xlim(-525, 0)

    plt.savefig('./figure-pdfs/f-dvdt-vs-vc.pdf')
    plt.show()


def plot_figure_dvdt_vs_vc():
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 5.5), sharex=True, sharey=True)
    fig.subplots_adjust(.1, .1, .95, .95)

    axs = axs.flatten()

    plot_dvdt_vs_vc(axs)

    for vals in [[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    plt.savefig('./figure-pdfs/f-dvdt-vs-vc.pdf')
    plt.show()


def plot_dvdt_vs_vc(ax, drug_name):
    meta = get_all_ap_meta()

    means_dvdt = []

    na_window = [1556.3, 1580]
    na_window = [2674.3, 2685]

    na_idx = np.array([int(val*25) for val in na_window])

    na_idx = np.array([int(val*25) for val in na_window])

    dat_dict = {'File': [],
               'Drug': [],
               'Baseline_AP': [],
               '3xEFPC_AP': [],
               '10xEFPC_AP': [],
               '20xEFPC_AP': [],
               'Washout_AP': [],
               'Baseline_VC': [],
               '3xEFPC_VC': [],
               '10xEFPC_VC': [],
               '20xEFPC_VC': [],
               'Washout_VC': []}

    for f, ap_meta in meta.items():
        if ap_meta['is_upstroke'].sum() < 45:
            continue

        # VC STUFF
        vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
        vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')

        vc_slice = vc_dat.iloc[na_idx[0]:na_idx[1], :]

        min_vals = vc_slice.min().values / vc_meta['cm'].values
        min_vals = [np.mean(min_vals[i*2:i*2+2]) for i in range(0, 5)]

        if np.min(min_vals) > -20:
            print(f'{f}: No sign of Na at any concentration')
            continue


        if 'flecainide' in f:
            dat_dict['Drug'].append('flecainide')
        else:
            dat_dict['Drug'].append('quinine')

        dat_dict['File'].append(f)

        # AP STUFF
        mean_dvdt = []
        for compound in ap_meta['compound'].unique():
            mean_dvdt.append([ap_meta[ap_meta['compound'] == compound]['dvdt_max'].mean(),
                               ap_meta[ap_meta['compound'] == compound]['dvdt_max'].std()])

        mean_dvdt = np.array(mean_dvdt)
        # USES A DIFFERENT CHANGE CALCULATION THAN f-flecainide_dvdt.plot_dvdt_changes(). 
        # IN THAT FUNCTION, I CALCULATE PRE- TO POST- WITHIN THE SAME TRIAL.
        # HERE, I LOOK AT DV/DT ACCROSS ALL TRIALS, REGARDLESS OF TIME.

        for i, compound in enumerate(['Baseline', '3xEFPC', '10xEFPC',
                                                  '20xEFPC', 'Washout']):
            dat_dict[f'{compound}_AP'].append(mean_dvdt[i,0])
            dat_dict[f'{compound}_VC'].append(min_vals[i])

    all_dat = pd.DataFrame(dat_dict)

    compounds = ['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Washout']
    is_changes = False 

    all_ap = []
    all_vc = []

    for i in range(0, 5):
        compound = compounds[i]

        # Current VC
        curr_vc_flec = all_dat[all_dat['File'].str.contains(drug_name)][f'{compounds[i]}_VC']
        # Next VC
        #next_vc_flec = all_dat[all_dat['File'].str.contains(drug_name)][f'{compounds[i+1]}_VC']

        # Current AP 
        curr_ap_flec = all_dat[all_dat['File'].str.contains(drug_name)][f'{compounds[i]}_AP']
        # Next AP 
        #next_ap_flec = all_dat[all_dat['File'].str.contains(drug_name)][f'{compounds[i+1]}_AP']

        if is_changes:
            # PLOT CHANGES
            ax.scatter(next_vc_flec-curr_vc_flec, next_ap_flec-curr_ap_flec, color='r', alpha=.4, label=f'{drug_name.capitalize()} ({compounds[i]})')
        else:
            # PLOT ABSOLUTES
            if compounds[i] == 'Baseline':
                ax.scatter(curr_vc_flec, curr_ap_flec, color=(i/4.0, 0, 0), alpha=.4, label=f'Baseline')
            elif 'EFPC' in compounds[i]:
                ax.scatter(curr_vc_flec, curr_ap_flec, color=(i/4.0, 0, 0), alpha=.4, label=f'{drug_name.capitalize()} ({compounds[i]})')
            else:
                ax.scatter(curr_vc_flec, curr_ap_flec, color='k', alpha=.4, marker='^', label=f'{drug_name.capitalize()} ({compounds[i]})')


        all_ap += curr_ap_flec.tolist()
        all_vc += curr_vc_flec.tolist()
        ax.legend(loc=3)

    slope, intercept, r_value, p_value, std_err = stats.linregress(
            all_vc, all_ap)
    print('VC vs AP')
    #print(f'\tFlat MDP: {np.mean(mdp_flat)} +/- {np.std(mdp_flat)}')
    print(f'p value of {p_value}')
    print(f'R value of {r_value}')


    regplot(all_vc,all_ap, color='grey', scatter=False, ax=ax, ci=None)
    from scipy.stats.stats import pearsonr
    feature_corrs = pearsonr(all_vc, all_ap)

    print(feature_corrs)


    import pdb
    pdb.set_trace()




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

    return [max_dvdt, is_na]


def get_all_ap_meta():
    #dmso_files = get_valid_cells('dmso')
    flecainide_files = get_valid_cells('flecainide')
    quinine_files = get_valid_cells('quinine')

    all_files = flecainide_files + quinine_files
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
    #plot_figure_dvdt_vs_vc()
    plot_single_dvdt_vs_vc()


if __name__ == '__main__':
    main()
