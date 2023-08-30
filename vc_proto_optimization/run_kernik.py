import myokit
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy.signal import find_peaks

from cell_models import kernik, protocols


def plot_kernik_NaL():
    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

    for i, mod_name in enumerate(
            ['mmt-files/kernik_2019_NaL.mmt', 'mmt-files/kernik_2019_mc.mmt']):
        mod, proto, x = myokit.load(mod_name)
        sim = myokit.Simulation(mod, proto)
        t = 3000 
        times = np.arange(0, t, .1)
        dat = sim.run(t, log_times=times)

        t = dat['engine.time']
        v = dat['membrane.V']

        if i == 0:
            axs[0].plot(t, v, label='With I_NaL')
            for i, curr in enumerate(['ina.i_Na', 'INaL.INaL']):
                axs[i+1].plot(t, dat[curr])
        else:
            axs[0].plot(t, v, label='No I_NaL')
            for i, curr in enumerate(['ina.i_Na']):
                axs[i+1].plot(t, dat[curr])


    fs = 18
    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    axs[0].set_ylabel('mV', fontsize=fs)
    axs[0].legend()
    axs[1].set_ylabel('I_Na', fontsize=fs)
    axs[2].set_ylabel('I_NaL', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)
    #ax.tick_params(labelsize=fs-4)
    #axs[0].set_xlim(0, 500)
    plt.show()


def compare_kernik_ord_sodium():
    proto = myokit.load_protocol('mmt-files/sodium_proto.mmt')
    proto_peak_v = list(range(-80, 60, 10))
    proto_peak_v = [-80, -70, -60, -50, -40, -35, -30, -25, -20, -15, -10, 0, 10, 20, 30, 40, 50] 

    mod = myokit.load_model('mmt-files/kernik_2019_mc.mmt')
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    t_max = proto.characteristic_time()

    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(t_max, log_times=times)

    kc_peaks = find_peaks(-np.array(dat['ina.i_Na']), distance=4000)
    ax.plot(proto_peak_v, np.array(dat['ina.i_Na'])[kc_peaks[0]], 'o-', label='Kernik')

    #axs[0].plot(dat['engine.time'], dat['membrane.V'])
    #axs[1].plot(dat['engine.time'], dat['ina.i_Na'], label='Kernik')

    mod = myokit.load_model('mmt-files/tor_ord_endo.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.v')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    t_max = proto.characteristic_time()

    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(t_max, log_times=times)

    tor_peaks = find_peaks(-np.array(dat['INa.INa']), distance=11000)

    ax.plot(proto_peak_v, np.array(dat['INa.INa'])[tor_peaks[0]], 'o-', label='Tor-ORd')

    #axs[1].plot(dat['engine.time'], dat['INa.INa'], label='Tor-ORd')
    #axs[1].set_xlim(995, 1015)

    ax.legend()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Voltage (mV)')
    ax.set_ylabel('INa current (A/F)')

    plt.show()


def compare_kernik_py_mmt_art():
    k_mod = kernik.KernikModel(is_exp_artefact=True)
    proto = protocols.VoltageClampProtocol([
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-80, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-70, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-60, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-50, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-40, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-35, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-30, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-25, 200),
                                    protocols.VoltageClampStep(-90, 2000),
                                    protocols.VoltageClampStep(-20, 200)
                                    ])
    tr = k_mod.generate_response(proto, is_no_ion_selective=False)

    mod = myokit.load_model('mmt-files/kernik_2019_NaL_art.mmt')
    proto = myokit.load_protocol('mmt-files/sodium_proto.mmt')
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    p = mod.get('engine.pace')
    p.set_binding(None)

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    t_max = proto.characteristic_time()

    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(t_max, log_times=times)

    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], np.array(dat['artifact.i_out'])/60, label='Kernik_mmt')

    axs[0].plot(tr.t, tr.y)
    axs[1].plot(tr.t, tr.current_response_info.get_current_summed(),
            label='Kernik_py')

    axs[0].set_ylabel('Membrane Voltage (mV)')
    axs[1].set_ylabel('I_out (pA/pF)')
    axs[1].set_xlabel('Time (ms)')

    axs[1].legend()
    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.show()


def compare_kernik_artifact():
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    mod = myokit.load_model('mmt-files/kernik_2019_NaL_art.mmt')
    proto = myokit.load_protocol('mmt-files/sodium_proto.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    t_max = proto.characteristic_time()

    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(t_max, log_times=times)

    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], np.array(dat['artifact.i_out'])/60, label='Kernik_mmt_NaL Artifact')


    mod = myokit.load_model('mmt-files/kernik_2019_NaL.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    t_max = proto.characteristic_time()

    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(t_max, log_times=times)

    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], label='Kernik_mmt_NaL')

    axs[0].set_ylabel('Membrane Voltage (mV)')
    axs[1].set_ylabel('I_out (pA/pF)')
    axs[1].set_xlabel('Time (ms)')

    axs[1].legend()
    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    plt.show()


#plot_tor_ord_baseline()
#plot_tor_ord_ikr()
#plot_tor_ord_ikr_ik1()
plot_kernik_NaL()
#compare_kernik_ord_sodium()
#compare_kernik_py_mmt_art()
#compare_kernik_artifact()
