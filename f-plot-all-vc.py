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


def plot_all_vc(channel=None):
    all_files = get_valid_cells('dmso') + get_valid_cells('flecainide') + get_valid_cells('quinine')
    #all_files = get_valid_cells('quinine')
    #all_files = get_valid_cells('dmso')
    #all_files = get_valid_cells('flecainide')
    

    windows = {'I_Na': [1550, 1580],
               'I_Kr': [3670, 3730],
               'I_to': [4856, 4876],
               'I_CaL': [4245, 4260]}
    ch_ylims = {'I_Na': None,
                'I_Kr': [-20, 8],
                'I_to': [-30, 10], 
                'I_CaL': [-30, 5]}

    if channel is not None:
        window = windows[channel]
        ylims = ch_ylims[channel]
    else:
        window = None
        window = [0, 10000]
        ylims = None

    fig, axs = plt.subplots(2, 1, figsize=(6.5, 5), sharex=True)
    fig.subplots_adjust(.12, .1, .95, .95)
    fig.suptitle(f'{channel} segment')#, fontsize=16)

    vc_pts = get_vc_pts()
    axs[0].plot(vc_pts[:, 0], vc_pts[:, 1], 'k')

    for i, f in enumerate(all_files):
        plot_window(f, window, axs[1], ylims=ylims)

    #plot_models(axs[0], axs[1], 'Kernik')
    #plot_models(axs[0], axs[1], 'Paci')

    axs[0].set_ylabel('Voltage (mV)')

    axs[1].set_ylabel('pA/pF')
    #axs[2].set_ylabel('dA_F/dt')
    axs[1].set_xlabel('Time (ms)')
    axs[0].set_xlim(window[0], window[1])
    #axs[1].legend()

    for ax in axs:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    if channel is None:
        plt.savefig(f'./figure-pdfs/f-all-vc.pdf')
    else:
        plt.savefig(f'./figure-pdfs/f-all-vc-{channel}.pdf')

    plt.show()


# PLOTTING FUNCTIONS
def plot_window(f, window, ax_c, ylims=None):
    vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
    vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')
    times = np.linspace(0, int(vc_dat.shape[0]/25), vc_dat.shape[0]+1)[0:-1]

    if window is not None:
        idx = np.array([int(val*25) for val in window])
        vc_slice = vc_dat.iloc[idx[0]:idx[1], :]
        times = times[idx[0]:idx[1]]
    else:
        vc_slice = vc_dat

    
    k = vc_dat.keys()[1]
    curr_dat = vc_slice[k].values
    curr_dat_1 = moving_average(curr_dat, n=5)
    times_1 = moving_average(times, n=5)

    curr_meta = vc_meta[vc_meta['sweep'] == k]
    ax_c.plot(times_1, curr_dat_1/vc_meta['cm'].mean(), label=f, color='grey', alpha=.4)


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
