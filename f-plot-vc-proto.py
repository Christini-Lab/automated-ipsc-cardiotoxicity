import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib

from utility import get_vc_pts

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)
matplotlib.rcParams['pdf.fonttype'] = 42



def plot_figure():
    fig, ax = plt.subplots(1, 1, figsize=(6, 3), sharex=True)
    fig.subplots_adjust(.13, .17, .95, .9)

    vc_pts = get_vc_pts()
    ax.plot(vc_pts[:, 0]-500, vc_pts[:, 1], 'k')
    ax.set_ylabel('Voltage (mV)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 9400)
    ax.set_xlabel('Time (ms)')

    label_dict = {'Na1': 520, 'Na2': 1057, 'NaL': 2225, 'Kr': 3197, 'CaL': 3749, 'to': 4365, 'K1': 4955, 'f': 6388, 'Ks': 8905}

    for name, time in label_dict.items():
        ax.axvline(time, color='grey', linestyle='--', alpha=.4)
        ax.text(time-50, 70, name, fontsize=9)

    plt.savefig(f'./figure-pdfs/f-vc-proto.pdf')

    plt.show()

plot_figure()
