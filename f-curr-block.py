import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib

from utility import get_valid_cells, get_vc_pts, get_model_response

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)
matplotlib.rcParams['pdf.fonttype'] = 42



def plot_drug_window2(f_name, window, concs, folder):
    fig, ax = plt.subplots(1, 1, figsize=(4, 3), sharex=True)
    fig.subplots_adjust(.2, .2, .95, .95)

    plot_window(f_name, window, ax, concs=concs)

    #vc_pts = get_vc_pts()
    #axs[0].plot(vc_pts[:, 0], vc_pts[:, 1], 'k')
    #axs[0].set_ylabel('Voltage (mV)')
    #axs[0].spines['top'].set_visible(False)
    #axs[0].spines['right'].set_visible(False)
    #axs[0].set_xlim(window[0], window[1])

    ax.set_ylabel('Current (A/F)')
    ax.set_xlabel('Time (ms)')
    ax.set_ylim(-120, 0)

    plt.savefig(f'./{folder}/f-{concs}.pdf')

    plt.show()


def plot_drug_window(drug, channel, concs):
    all_files = get_valid_cells(drug) 
    windows = {'I_Na': [1550, 1580],
               'I_Kr': [3670, 3730],
               'I_to': [4856, 4876],
               'I_CaL': [4245, 4260]}
    ch_ylims = {'I_Na': None,
                'I_Kr': [-20, 8],
                'I_to': [-30, 10], 
                'I_CaL': [-30, 5]}

    window = windows[channel]
    ylims = ch_ylims[channel]

    fig, axs = plt.subplots(len(all_files)+1, 1, figsize=(6.5, 5), sharex=True)
    fig.subplots_adjust(.12, .1, .95, .95)
    fig.suptitle(f'{drug} Effect on {channel} segment')#, fontsize=16)
    
    for i, f in enumerate(all_files):
        if i == 0:
            plot_window(f, window, axs[i+1], axs[0], ylims=ylims)
        else:
            plot_window(f, window, axs[i+1], ylims=ylims)

    axs[-1].set_xlabel('Time (ms)')

    plt.savefig(f'./figure-pdfs/f-{drug}-{channel}.pdf')
    plt.show()


# PLOTTING FUNCTIONS
def plot_window(f, window, ax_c, concs):
    idx = np.array([int(val*25) for val in window])

    vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
    vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')

    vc_slice = vc_dat.iloc[idx[0]:idx[1], :]

    times = np.linspace(0, int(vc_dat.shape[0]/25), vc_dat.shape[0]+1)[0:-1]
    times = times[idx[0]:idx[1]]

    mask = [v[2]=='2' for v in vc_meta['sweep'].str.split('_')]

    for i, k in enumerate(vc_dat.keys()):
        if mask[i]:
            continue
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

        ax_c.plot(times, curr_dat/vc_meta['cm'].mean(), c=cols)

        if concs == curr_meta['compound'].values[0]:
            break

    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)
    ax_c.set_ylabel('A/F')
    ax_c.set_xlim(times[0], times[-1])

    #scales = {'Baseline': {'g_Kr': 1, 'g_to': 1, 'g_Na': 1, 'g_CaL': 1},
    #          '3xEFPC': {'g_Kr': .3, 'g_to': .75, 'g_Na': .9, 'g_CaL': .98},
    #          '10xEFPC': {'g_Kr': .85, 'g_to': .45, 'g_Na': .52, 'g_CaL': .87},
    #          '20xEFPC': {'g_Kr': .9, 'g_to': .6, 'g_Na': .8, 'g_CaL': .3}
    #          }

    #for conc, scale_vals in scales.items():
    #    dat = get_model_response('Kernik', scale_vals)


def main():
    window = [1555, 1560]
    file = 'flecainide_221028_008_4'
    for conc in ['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Washout']:
        plot_drug_window2(file, window, concs=conc, folder='figure-pdfs/f-increased-flec-conc')

    #window = [3675, 3730]
    #file = 'flecainide_221028_003_4'
    ##file = 'flecainide_220901_009_3'
    #for conc in ['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Washout']:
    #    plot_drug_window2(file, window, concs=conc, folder='figure-pdfs/f-increased-flec-conc-ikr')



if __name__ == '__main__':
    main()


