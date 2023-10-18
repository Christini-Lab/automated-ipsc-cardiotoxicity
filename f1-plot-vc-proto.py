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
    ax.plot(vc_pts[:, 0], vc_pts[:, 1], 'k')
    ax.set_ylabel('Voltage (mV)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 9800)
    ax.set_xlabel('Time (ms)')

    #label_dict = {'Na1': 520, 'Na2': 2180, 'NaL': 2225, 'Kr': 3197, 'CaL': 3749, 'to': 4365, 'K1': 4955, 'f': 6388, 'Ks': 8905}
    label_dict = {'Na': 2180, 'NaL': 2225, 'Kr': 3197, 'CaL': 3749, 'to': 4365, 'K1': 4955, 'f': 6380, 'Ks': 8900}
    for name, time in label_dict.items():
        time += 500
        ax.axvline(time, color='grey', linestyle='--', alpha=.4)
        if name == 'Na':
            ax.text(time-350, 60, name, fontsize=9)    
        elif name == 'NaL':
            ax.text(time+50, 60, name, fontsize=9)  
        else: 
            ax.text(time-50, 60, name, fontsize=9)

    plt.savefig(f'./figure-pdfs/f-vc-proto.pdf')

    plt.show()

plot_figure()
