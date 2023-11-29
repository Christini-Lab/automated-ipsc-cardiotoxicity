import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utility import get_valid_cells, get_vc_pts
from utility_funcs import get_cell_objects

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)


def plot_figure():
    fig = plt.figure(figsize=(6.5, 3))
    fig.subplots_adjust(.1, .2, .95, .95)

    grid = fig.add_gridspec(1, 3, hspace=.4, wspace=.4)
                                #height_ratios=[1,2])

    axs = []

    sub = grid[2]
    corr_sub = sub.subgridspec(1, 1)
    ax = fig.add_subplot(corr_sub[0])
    axs.append(ax)
                            
    plot_ikr_segment(axs[0])
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Current (A/F)')
    ax.axvline(3696, linestyle='--', color='grey')

    sub = grid[:, 1]
    corr_sub = sub.subgridspec(2, 1, height_ratios=[1,2])
    ax_v = fig.add_subplot(corr_sub[0])
    ax_i = fig.add_subplot(corr_sub[1])
    #8, 13 
    window = [3680, 3710]
    cell_name = get_valid_cells('quinine')[8]
    plot_apc_drug_ex(ax_i, cell_name, window)
    ax_i.text(3688, 5, 'Automated')
    ax_v.plot([3680, 3684, 3684, 3698, 3698, 3710], 
                [24, 24, -37, -37, 28, 28], color='k')
    ax_v.set_ylim(-40, 60)
    ax_v.set_axis_off()
    #ax_v.text("24 mV", xy=[3680, 24])
    ax_v.text(3680, 28, "24 mV")
    ax_v.text(3685, -33, "-37 mV")
    ax_v.text(3699, 32, "28 mV")
    axs.append(ax_i)
    axs.append(ax_v)
    ax_i.set_xlabel('Time (ms)')
    ax_i.set_ylabel('Current (A/F)')
    ax_i.annotate("", xy=(3696, -3), xytext=(3696, -8), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='k'))


    sub = grid[:, 0]
    corr_sub = sub.subgridspec(2, 1, height_ratios=[1,2])
    ax_v = fig.add_subplot(corr_sub[0])
    ax_i = fig.add_subplot(corr_sub[1])
    plot_manual_ex(ax_i)
    ax_v.plot([855, 857, 857, 865, 865, 870], 
                [6, 6, -41, -41, 8.5, 8.5], color='k')
    ax_v.set_ylim(-42, 20)
    ax_v.set_axis_off()
    #ax_v.text("24 mV", xy=[3680, 24])
    ax_v.text(855, 9, "6 mV")
    ax_v.text(857.5, -38, "-41 mV")
    ax_v.text(866, 11, "8 mV")
    axs.append(ax_i)
    axs.append(ax_v)
    ax_i.set_xlabel('Time (ms)')
    ax_i.annotate("", xy=(863, -1.5), xytext=(863, -6), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='k'))

    for i, ax in enumerate(axs):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    axs[0].set_title('C', y=.94, x=-.35)
    axs[2].set_title('B', y=.9, x=-.15)
    axs[4].set_title('A', y=.9, x=-.15)


    plt.savefig('figure-pdfs/f-no-ikr.pdf')

    plt.show()


def plot_ikr_segment(ax, curr_name='I_Kr'):
    window = [3680, 3710]
    ax.set_xlim(window[0], window[1])
    plot_i(ax, window, curr_name='I_Kr')
    ax.set_ylim(-25, 10)


def plot_apc_drug_ex(ax, cell_name, window):
    labs = ['Baseline', None, '3xEFPC', None, '10xEFPC', None, '20xEFPC', None, 'Wash', None]
    vc_dat = pd.read_csv(f'data/cells/{cell_name}/vc_df.csv')
    vc_meta = pd.read_csv(f'data/cells/{cell_name}/vc_meta.csv')

    times = np.linspace(0, int(vc_dat.shape[0]/25), vc_dat.shape[0]+1)[0:-1]
    idx = np.array([int(val*25) for val in window])

    i=0
    for k in ['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Washout']:
        curr_meta = vc_meta[vc_meta['compound'] == k]

        curr_dat = vc_dat[curr_meta['sweep'].values[0]].values
        #curr_dat_2 = vc_dat[curr_meta['sweep'].values[1]].values

        mv_avg = 3 

        #curr_dat = np.mean([curr_dat_1, curr_dat_2], 0)
        curr_slice = curr_dat[idx[0]:idx[1]]
        curr_dat = moving_average(curr_slice, n=mv_avg)/vc_meta['cm'].mean()

        if 'Baseline' == k:
            cols = c=(0, 0, 0)
        if '3xEFPC' == k:
            cols = c=(1, .7, .7)
        if '10xEFPC' == k:
            cols =(1, .4, .4)
        if '20xEFPC' == k:
            cols =(1, 0, 0)
        if 'Washout' == k:
            cols =(.5, .5, .5)

        it = window
        st = int(it[0]/times[1])
        end = int(it[1]/times[1])
        t = times[st:end]
        t = moving_average(t, n=mv_avg)
        #curr = curr_dat[st:end]#*1E9
        ax.plot(t, curr_dat, c=cols, label=labs[i])
        i+=1
    

def plot_manual_ex(ax_i):
    pre_drug = pd.read_csv('./data/manual_data/pre-drug_vc_proto.csv')
    post_drug = pd.read_csv('./data/manual_data/post-drug_vc_proto.csv') 

    start, end = 8550, 8700

    #fig, axs = plt.subplots(2, 1)
    #axs[0].plot(pre_drug['Time (s)'], pre_drug['Voltage (V)'])
    t = pre_drug['Time (s)'][start:end] *1000
    ax_i.plot(t, pre_drug['Current (pA/pF)'][start:end], 'k')
    ax_i.plot(t, post_drug['Current (pA/pF)'][start:end], 'r')
    ax_i.set_ylim(-8, 8)
    ax_i.text(855, 5.5, 'Manual (Clark, 2022)')
    #ax_v.plot(t, post_drug['Voltage (V)'][start:end]*1000, color='k')

    #ax.plo




def plot_i(ax, window=[0, 9444], curr_name=None):
    all_files = get_valid_cells('dmso') + get_valid_cells('flecainide') + get_valid_cells('quinine')
    all_files = all_files[2:22]

    idx = np.array([int(val*25) for val in window])

    all_dat = []

    mv_avg = 5
    if window[1] - window[0] < 20:
        mv_avg = 2

    for i, f in enumerate(all_files):
        vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
        vc_slice = vc_dat.iloc[idx[0]:idx[1], :]
        vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')
        k = vc_dat.keys()[1]
        curr_dat = vc_slice[k].values
        curr_dat_1 = moving_average(curr_dat, n=mv_avg)/vc_meta['cm'].mean()

        all_dat.append(curr_dat_1)
        print(i)


    times = np.linspace(0, int(vc_dat.shape[0]/25), vc_dat.shape[0]+1)[0:-1]
    times = times[idx[0]:idx[1]]
    times_1 = moving_average(times, n=mv_avg)

    arr_indices = [0, 1, 2]
    arr_indices = list(range(0, len(all_dat)))

    for idx in arr_indices:
        ax.plot(times_1, all_dat[idx], color='grey', alpha=.3)


    if curr_name == 'I_Na':
        min_idx = np.argmin(np.min(all_dat, 1))
        max_idx = np.argmax(np.min(all_dat, 1))

        ax.plot(times_1, all_dat[min_idx])
        ax.plot(times_1, all_dat[max_idx])

    if curr_name == 'I_f':
        min_idx = np.argmin(np.min(all_dat, 1))
        max_idx = np.argmax(np.min(all_dat, 1))

        idxs = [2, 3]
        #idxs = [2, 4, 5, 6]

        [ax.plot(times_1, all_dat[idx]) for idx in idxs]

        #ax.plot(times_1, all_dat[min_idx])
        #ax.plot(times_1, all_dat[max_idx])
    
    if curr_name == 'I_Kr':
        min_idx = np.argmin(np.min(all_dat, 1))
        max_idx = np.argmax(np.min(all_dat, 1))

        idxs = [2, 3]
        #idxs = [2, 4, 5, 6]

        [ax.plot(times_1, all_dat[idx]) for idx in idxs]

        #ax.plot(times_1, all_dat[min_idx])
        #ax.plot(times_1, all_dat[max_idx])


# PLOTTING FUNCTIONS
def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]

    return np.array(new_vals)



def main():
    plot_figure()

if __name__ == "__main__":
    main()
