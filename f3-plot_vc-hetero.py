import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import myokit

from utility import get_valid_cells, get_vc_pts

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)


def plot_all_vc():
    fig = plt.figure(figsize=(6.5, 6))
    fig.subplots_adjust(.1, .12, .95, .95)

    axs = []
    axs_v = []

    #PANEL A
    grid = fig.add_gridspec(2, 2, hspace=.3, wspace=.2)#, height_ratios=[1,2,1,1])
    sub = grid[0, 0:]
    corr_subgrid = sub.subgridspec(2, 1, hspace=.1)
    ax_v = fig.add_subplot(corr_subgrid[0])
    ax_i = fig.add_subplot(corr_subgrid[1])

    window = [0, 9444]
    vc_pts = get_vc_pts()
    ax_v.plot(vc_pts[:, 0], vc_pts[:, 1], 'k')
    ax_v.set_xlim(window[0], window[1])
    ax_i.set_xlim(window[0], window[1])
    plot_i(ax_i, window)
    ax_i.set_ylim(-40, 40)

    ax_i.set_xlabel('Time (ms)')
    ax_v.set_ylabel('mV')
    ax_i.set_ylabel(r'$I_{out}$ (A/F)')
    ax_v.tick_params(labelbottom=False)

    label_dict = {'Na': 2180, 'NaL': 2225, 'Kr': 3197, 'CaL': 3749, 'to': 4365, 'K1': 4955, 'f': 6380, 'Ks': 8900}
    ax = ax_v
    for name, time in label_dict.items():
        time += 500
        ax.axvline(time, color='grey', linestyle='--', alpha=.4)
        if name == 'Na':
            ax.text(time-350, 60, name, fontsize=9)    
        elif name == 'NaL':
            ax.text(time+50, 60, name, fontsize=9)  
        else: 
            ax.text(time-50, 60, name, fontsize=9)
    

    axs.append(ax_v)
    axs.append(ax_i)
    axs_v.append(ax_v)

    windows = [[2670, 2685], [5600, 7000]]
    ylims_v = [[-120, -10], [-125, -50]]
    ylims_i = [[-510, 0], [-30, 0]]

    curr_names = ['I_Na', 'I_f']
    panel_titles = [r'$I_{Na}$', r'$I_{f}$']
    x_ticks = [[2670, 2675, 2680, 2685], 
                    [5600, 6000, 6400, 6800]]


    #PANELS B-C
    for i in [0, 1]:
        window = windows[i]

        sub = grid[1, i]
        corr_subgrid = sub.subgridspec(2, 1, hspace=.1)
        ax_v = fig.add_subplot(corr_subgrid[0])
        ax_i = fig.add_subplot(corr_subgrid[1])

        ax_v.plot(vc_pts[:, 0], vc_pts[:, 1], 'k')
        ax_v.set_xlim(window[0], window[1])
        ax_i.set_xlim(window[0], window[1])

        plot_i(ax_i, window, curr_name=curr_names[i])

        ax_i.set_xlabel('Time (ms)')
        if i == 0:
            ax_v.set_ylabel('mV')
            ax_i.set_ylabel(r'$I_{out}$ (A/F)')

        ax_v.tick_params(labelbottom=False)
        ax_v.set_ylim(ylims_v[i][0], ylims_v[i][1])
        ax_i.set_ylim(ylims_i[i][0], ylims_i[i][1])

        ax_i.set_xticks(x_ticks[i])

        ax_v.set_title(panel_titles[i], y=.8)

        axs.append(ax_v)
        axs.append(ax_i)
        axs_v.append(ax_v)

    letter_label = ['A', 'B', 'C']
    letter_coords = [[-700, 50], [2669, -10], [5400, -50]]

    for i, ax in enumerate(axs):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for i, ax in enumerate(axs_v):
        ax.text(letter_coords[i][0], letter_coords[i][1], letter_label[i], fontsize=11)


    plt.savefig('figure-pdfs/f-vc_hetero.pdf')
    plt.show()


def plot_i(ax, window=[0, 9444], curr_name=None):
    all_files = get_valid_cells('dmso') + get_valid_cells('flecainide') + get_valid_cells('quinine')
    all_files = all_files[0:25]

    idx = np.array([int(val*25) for val in window])

    all_dat = []

    mv_avg = 5
    if window[1] - window[0] < 20:
        mv_avg = 2

    all_cap = []
    all_rseries = []
    all_seal = []

    for i, f in enumerate(all_files):
        vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
        vc_slice = vc_dat.iloc[idx[0]:idx[1], :]
        vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')
        all_cap.append(vc_meta['cm'].mean())
        all_rseries.append(vc_meta['r_series'].mean())
        all_seal.append(vc_meta['r_seal'].mean())
        k = vc_dat.keys()[1]
        curr_dat = vc_slice[k].values
        curr_dat_1 = moving_average(curr_dat, n=mv_avg)/vc_meta['cm'].mean()

        all_dat.append(curr_dat_1)
        print(i)


    #Measure capacitance AND rseries
    np.mean(all_cap)/1E-12
    np.std(all_cap)/1E-12
    np.mean(all_rseries)/1E6
    np.std(all_rseries)/1E6

    import pdb
    pdb.set_trace()

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
 

# PLOTTING FUNCTIONS
def plot_models(ax_v, ax_i, model_name):
    if model_name == 'Kernik':
        mod = myokit.load_model(
                './vc_proto_optimization/mmt-files/kernik_2019_NaL_art.mmt')
        styles = 'r--'

        cm = 20

        mod['geom']['Cm'].set_rhs(cm)
        mod['artifact']['r_series'].set_rhs(10E-3)

        p = mod.get('engine.pace')
        p.set_binding(None)

        v_cmd = mod.get('artifact.v_cmd')
        v_cmd.set_rhs(0)
        v_cmd.set_binding('pace') # Bind to the pacing mechanism

        proto = myokit.Protocol()

        vc_steps = get_vc_pts()

        for i in range(0, int(len(vc_steps)/2)):
            curr_v = vc_steps[i*2][1]
            curr_duration = vc_steps[2*i+1][0] - vc_steps[2*i][0]
            
            proto.add_step(curr_v, curr_duration)


        t_max = proto.characteristic_time()

        times = np.arange(0, t_max, 0.1)

        sim = myokit.Simulation(mod, proto)
        dat = sim.run(t_max, log_times=times)

        ax_i.plot(dat['engine.time'], np.array(dat['artifact.i_out']) / cm, 'r--', label='Kernik_mmt_NaL Artifact', linewidth=1.5)
        ax_v.plot(dat['engine.time'], dat['membrane.V'], 'r--')

    else:
        mod = myokit.load_model(
                './vc_proto_optimization/mmt-files/paci_cardio_lei.mmt')
        styles = 'b--'

        mod['voltageclamp']['rseries'].set_rhs(10E6)
        mod['voltageclamp']['rseries_est'].set_rhs(10E6)

        proto = myokit.Protocol()
        vc_steps = get_vc_pts()

        for i in range(0, int(len(vc_steps)/2)):
            curr_v = vc_steps[i*2][1]
            curr_duration = vc_steps[2*i+1][0] - vc_steps[2*i][0]
            
            proto.add_step(curr_v/1000, curr_duration/1000)

        t_max = proto.characteristic_time()

        times = np.arange(0, t_max, 0.1/1000)

        sim = myokit.Simulation(mod, proto)
        dat = sim.run(t_max, log_times=times)

        ax_i.plot(times*1000, np.array(dat['voltageclamp.Iout'])*1E12, 'b--', label='Paci Artifact', linewidth=1.5)
        ax_v.plot(times*1000, np.array(dat['Membrane.Vm'])*1000, 'b--')


def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]

    return np.array(new_vals)


def main():
    plot_all_vc()


if __name__ == '__main__':
    main()
