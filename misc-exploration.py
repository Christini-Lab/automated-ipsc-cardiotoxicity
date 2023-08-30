import pandas as pd
import numpy as np
from os import listdir
import matplotlib.pyplot as plt



def plot_vc_dat(cell_folder):
    labs = ['Baseline', None, '3xEFPC', None, '10xEFPC', None, '20xEFPC', None, 'Wash', None]

    pts = get_vc_pts()

    vc_dat = pd.read_csv(f'data/cells/{cell_folder}/vc_df.csv')
    vc_meta = pd.read_csv(f'data/cells/{cell_folder}/vc_meta.csv')


    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    i = 0
    for k in vc_dat.keys():
        times = np.linspace(0, int(vc_dat.shape[0]/25), vc_dat.shape[0]+1)[0:-1]
        curr_dat = vc_dat[k].values

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

        cm = vc_meta['cm'].mean()
    
        curr_dat = moving_average(curr_dat, n=10)
        times = moving_average(times, n=10)


        axs[1].plot(times, curr_dat/cm, c=cols, label=labs[i])

    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)


    pts = np.array(pts)

    axs[0].plot(pts[:, 0], pts[:, 1], 'k')

    axs[1].set_xlabel('Time (ms)')
    axs[0].set_ylabel('Voltage (mV)')
    axs[1].set_ylabel('Current (A/F)')

    axs[1].legend()
    plt.show()


def get_vc_pts():
    steps = {
            'Voltage':  [-105, -36, -80, -112, -2, -116, -31, -80, 24, -37, 28, -80, 3, -90, -80, -68, 57, -80, -68, -120, -80, -120, -120, -77, -80, 23, 40, -1, -80, -50, 0],
            'Duration': [12, 40, 500, 5, 364, 754, 88, 500, 422, 14, 40, 500, 10, 40, 500, 71, 43, 500, 45, 46, 500, 18, 876, 40, 500, 357, 1620, 40, 500, 600, 100]
            }

    pts = [[0, -80]]
    t = 999 
    curr_voltage = -80

    for i in range(0, len(steps['Voltage'])):
        pts.append([t, curr_voltage])
        curr_voltage = int(float(steps['Voltage'][i]))
        pts.append([t, curr_voltage])

        t += steps['Duration'][i] 

    return np.array(pts)


def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]
    return np.array(new_vals)


plot_vc_dat('quinine_221206_009_1')
