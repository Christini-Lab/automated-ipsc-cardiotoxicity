import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import pointplot, swarmplot, histplot
import pickle
from scipy.stats import ttest_ind
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



def plot_figure():
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6))
    fig.subplots_adjust(.1, .15, .95, .95, wspace=.3)

    axs = axs.flatten()

    plot_aps(axs[0])

    plot_ap_features(axs[1], 'APD90')
    plot_ap_features(axs[2], 'MDP')
    plot_ap_features(axs[3], 'dVdt')

    axs[1].set_xlabel(r'$APD_{90}$ (ms)')
    axs[2].set_xlabel(r'MDP (mV)')
    axs[3].set_xlabel(r'$dV/dt_{max}$ (V/s)')

    for vals in [[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    plt.savefig(f'./figure-pdfs/f-ap_hetero.pdf')
    plt.show()


def plot_aps(ax):
    all_files = get_valid_cells('dmso') + get_valid_cells('flecainide') + get_valid_cells('quinine')

    times = np.linspace(0, 1000, 25025) - 50

    st_idxs = 1000
    num_idxs = 4700 

    print(f'There are {len(all_files)} files.')

    for f in all_files:
        ap_dat = pd.read_csv(f'./data/cells/{f}/ap_df.csv')

        ax.plot(times[st_idxs:num_idxs], ap_dat.iloc[:, 4].values[st_idxs:num_idxs], color='k', alpha=.3)

    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Voltage (mV)')
    

def plot_ap_features(ax, feature):
    all_files = get_valid_cells('dmso') + get_valid_cells('flecainide') + get_valid_cells('quinine')

    times = np.linspace(0, 1000, 25025) - 50
    num_idxs = 5000

    print(f'There are {len(all_files)} files.')

    all_vals = []

    for f in all_files:
        ap_dat = pd.read_csv(f'./data/cells/{f}/ap_df.csv')

        if feature == 'APD90':
            feature_val = np.mean([get_apd(ap_dat.values[:, i])[0] for i in range(0, 5)]) 
        elif feature == 'dVdt':
            feature_val = np.mean([get_dvdt(ap_dat.values[:, i])[0] for i in range(0, 5)]) 
        else:
            feature_val = np.mean([get_mdp(ap_dat.values[:, i]) for i in range(0, 5)]) 

        all_vals.append(feature_val)


    print(f'The mean and SD of {feature} is {np.mean(all_vals)} +- {np.std(all_vals)}')
    histplot(all_vals, ax=ax, color='k')

    #ax.set_xlabel('Time (ms)')
    #ax.set_ylabel('Voltage (mV)')
    

def get_apd(ap):
    stim_region = ap[int(48*25):int(65*25)]

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


def get_dvdt(ap):
    stim_region = ap[int(48*25):int(65*25)]

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


def get_mdp(ap):
    return np.mean(ap[0:50])


def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]
    return np.array(new_vals)


plot_figure()
