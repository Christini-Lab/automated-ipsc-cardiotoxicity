import pandas as pd
import numpy as np
from os import listdir
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rc('legend', fontsize = 8)
matplotlib.rcParams['pdf.fonttype'] = 42



def save_all_vc_ap_data():
    all_cells = listdir(f'./data/cells')

    all_cells = [c for c in all_cells if 'DS_' not in c]

    for cell_folder in all_cells:
        if 'archive' in cell_folder:
            continue
        if '230412' not in cell_folder:
            continue
        #if 'flecainide_221028' not in cell_folder:
        #    continue
        print(cell_folder)
        save_ap_dat(cell_folder)

        save_vc_dat(cell_folder, None)
        #save_vc_dat(cell_folder, ['INa1', [1010, 1020]])
        #save_vc_dat(cell_folder, ['INa2', [1555, 1565]])
        #save_vc_dat(cell_folder, ['IKs', [9000, 9020]])


def save_ap_dat(cell_folder):
    ap_dat = pd.read_csv(f'data/cells/{cell_folder}/ap_df.csv')
    ap_meta = pd.read_csv(f'data/cells/{cell_folder}/ap_meta.csv')

    trial_nums = np.sort(np.unique([int(trial.split('_')[1]) for trial in ap_dat.keys()]))

    fig, axs = plt.subplots(1, 4, figsize=(16, 5), sharey=True, sharex=True)

    times = np.linspace(-.05, .951, ap_dat.shape[0]+1)[0:-1]

    labs = ['Baseline', None, '3xEFPC', None, '10xEFPC', None, '20xEFPC', None, 'Wash', None]

    for k in ap_dat.keys():
        curr_dat = ap_dat[k].values
        ax_index = np.where(int(k.split('_')[1]) == trial_nums)[0][0]
        curr_compound = ap_meta[ap_meta['sweep'] == k]

        curr_meta = ap_meta[ap_meta['sweep'] == k]

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

        axs[ax_index].plot(times[0:10000]*1000, curr_dat[0:10000], c=cols)


    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_xlabel('Time (ms)')

    axs[0].set_title('Baseline to 3xEFPC')
    axs[1].set_title('3xEFPC to 10xEFPC')
    axs[2].set_title('10xEFPC to 20xEFPC')
    axs[3].set_title('20xEFPC to DMSO wash')
    axs[0].set_ylabel('Voltage (mV)')

    axs[0].set_xlim(-10, 150)

    plt.savefig(f'data/cells/{cell_folder}/ap_dat.pdf')

    axs[0].set_xlim(-3, 7)

    plt.savefig(f'data/cells/{cell_folder}/ap_upstroke.pdf')

    plt.close()


def save_vc_dat(cell_folder, time_window, with_mv_avg=False):
    labs = ['Baseline', None, '3xEFPC', None, '10xEFPC', None, '20xEFPC', None, 'Wash', None]

    vc_file = open('./vc_proto_optimization/results/exp_14/long_proto_steps.csv', 'r')
    steps = vc_file.readlines()

    pts = [[0, -80]]
    t = 999
    curr_voltage = -80
    new_times = [999, 1011, 1051, 1551, 1556]

    for i, step in enumerate(steps):
        vals = step.split()

        if i < len(new_times):
            pts.append([new_times[i], curr_voltage])
            curr_voltage = int(float(vals[3][:-1]))
            pts.append([new_times[i], curr_voltage])
        else:
            pts.append([t, curr_voltage])
            curr_voltage = int(float(vals[3][:-1]))
            pts.append([t, curr_voltage])

        #print(f'Step {i}: Voltage of {}, and time of{}') 
        t += int(float(vals[-1]))

    vc_dat = pd.read_csv(f'data/cells/{cell_folder}/vc_df.csv')
    vc_meta = pd.read_csv(f'data/cells/{cell_folder}/vc_meta.csv')

    times = np.linspace(0, int(vc_dat.shape[0]/25), vc_dat.shape[0]+1)[0:-1]

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    i = 0
    for k in ['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Washout']:
    #for k in vc_dat.keys():
        curr_meta = vc_meta[vc_meta['compound'] == k]

        try:
            curr_dat_1 = vc_dat[curr_meta['sweep'].values[0]].values
        except:
            import pdb
            pdb.set_trace()
        curr_dat_2 = vc_dat[curr_meta['sweep'].values[1]].values

        curr_dat = np.mean([curr_dat_1, curr_dat_2], 0)

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

        if time_window is None:
            cm = curr_meta['cm'].mean()
            if with_mv_avg:
                times_new = moving_average(times, 10)
                curr_new = moving_average(curr_dat, 10)
                axs[1].plot(times_new, curr_new/cm, c=cols, label=k)
            else:
                axs[1].plot(times, curr_dat/cm, c=cols, label=labs[i])
            #axs[1].plot(times, curr_dat, c=cols, label=labs[i])
        else:
            it = time_window
            st = int(it[1][0]/times[1])
            end = int(it[1][1]/times[1])
            t = times[st:end]
            curr = curr_dat[st:end]*1E9
            axs[1].plot(t, curr, c=cols, label=labs[i])
        i+=1

    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)


    pts = np.array(pts)

    axs[0].plot(pts[:, 0], pts[:, 1], 'k')

    if time_window is not None:
        axs[0].set_xlim(time_window[1][0], time_window[1][1])
        axs[0].set_title(time_window[1])

    axs[1].set_xlabel('Time (ms)')
    axs[0].set_ylabel('Voltage (mV)')
    axs[1].set_ylabel('Current (A/F)')

    axs[1].legend()
    if time_window is None:
        import pdb
        pdb.set_trace()
        plt.savefig(f'data/cells/{cell_folder}/vc_dat.pdf',)
    else:
        plt.savefig(f'data/cells/{cell_folder}/{time_window[0]}_vc_dat.pdf',)

    curr_windows = {'INa1': [1000, 1050],
                    'INa2': [1550, 1580],
                    'INaL': [2650, 2800],
                    'IKr': [3660, 3750],
                    'ICaL': [4825, 5000],
                    'IK1': [5900, 6020],
                    'Ito': [4825, 5000],
                    'If': [5000, 8000],
                    'IKs': [8000, 9400]}

    for curr, window in curr_windows.items():
        axs[0].set_xlim(window[0], window[1])
        fig.suptitle(curr)
        st, end = int(window[0]*25), int(window[1]*25) 
        mins = [np.min(l.get_data()[1][st:end]) for l in axs[1].get_lines()]
        maxs = [np.max(l.get_data()[1][st:end]) for l in axs[1].get_lines()]

        min_val = np.min(mins)
        max_val = np.max(maxs)
        if 'Na' in curr:
            axs[1].set_ylim(min_val, 30)
        else:
            if min_val < -100:
                min_val = -100
            if max_val > 30:
                max_val = 30 

            axs[1].set_ylim(min_val, max_val)

        #plt.savefig(f'data/cells/{cell_folder}/{curr}_vc_dat.pdf')


    plt.close()


def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]

    return np.array(new_vals)



#save_all_vc_ap_data()

save_vc_dat('quinine_221206_009_1', None, with_mv_avg=True)
