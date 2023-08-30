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


def plot_figure_dvdt_vs_vc():
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 5.5), sharex=True, sharey=True)
    fig.subplots_adjust(.1, .1, .95, .95)

    axs = axs.flatten()

    plot_dvdt_vs_vc(axs)

    for vals in [[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    plt.savefig('./figure-pdfs/f-quinine-dvdt-vs-vc.pdf')
    plt.show()


def plot_dvdt_vs_vc(axs):
    meta = get_all_ap_meta()

    flecainide_means_dvdt = []

    na_window = [1556.3, 1580]
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


        if 'dmso' in f:
            dat_dict['Drug'].append('DMSO')
        else:
            dat_dict['Drug'].append('Flecainide')

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

    for i, ax in enumerate(axs):
        compound = compounds[i]

        # Current VC
        curr_vc_dmso = all_dat[all_dat['File'].str.contains('dmso')][f'{compounds[i]}_VC']
        curr_vc_flec = all_dat[all_dat['File'].str.contains('quinine')][f'{compounds[i]}_VC']
        # Next VC
        next_vc_dmso = all_dat[all_dat['File'].str.contains('dmso')][f'{compounds[i+1]}_VC']
        next_vc_flec = all_dat[all_dat['File'].str.contains('quinine')][f'{compounds[i+1]}_VC']

        # Current AP 
        curr_ap_dmso = all_dat[all_dat['File'].str.contains('dmso')][f'{compounds[i]}_AP']
        curr_ap_flec = all_dat[all_dat['File'].str.contains('quinine')][f'{compounds[i]}_AP']
        # Next AP 
        next_ap_dmso = all_dat[all_dat['File'].str.contains('dmso')][f'{compounds[i+1]}_AP']
        next_ap_flec = all_dat[all_dat['File'].str.contains('quinine')][f'{compounds[i+1]}_AP']

        if is_changes:
            # PLOT CHANGES
            ax.scatter(next_vc_dmso-curr_vc_dmso, next_ap_dmso-curr_ap_dmso, color='k', alpha=.4,
                    label=f'DMSO ({compounds[i]})')
            ax.scatter(next_vc_flec-curr_vc_flec, next_ap_flec-curr_ap_flec, color='r', alpha=.4,
                    label=f'Flecainide ({compounds[i]})')
        else:
            # PLOT ABSOLUTES
            ax.scatter(curr_vc_dmso, curr_ap_dmso, color='k', alpha=.4, label=f'DMSO ({compounds[i]})')
            ax.scatter(curr_vc_flec, curr_ap_flec, color='r', alpha=.4, label=f'Flecainide ({compounds[i]})')
        ax.legend(loc=3)

    if is_changes:
        axs[0].set_ylabel(r'$\Delta$$dV/dt_{max}$ (V/s)')
        axs[2].set_ylabel(r'$\Delta$$dV/dt_{max}$ (V/s)')

        axs[2].set_xlabel(r'$\Delta$Current (A/F)')
        axs[3].set_xlabel(r'$\Delta$Current (A/F)')
    else:
        axs[0].set_ylabel(r'$dV/dt_{max}$ (V/s)')
        axs[2].set_ylabel(r'$dV/dt_{max}$ (V/s)')

        axs[2].set_xlabel(r'Current (A/F)')
        axs[3].set_xlabel(r'Current (A/F)')



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
    plot_figure_dvdt_vs_vc()


if __name__ == '__main__':
    main()
