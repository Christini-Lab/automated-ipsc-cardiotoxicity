import pickle
from vc_opt_ga import simulate_model
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from utility_classes import VCProtocol, VCSegment
from os import listdir
import myokit



def get_seal_comparison(hold_val, is_shown=True):
    cm = 60

    mod = myokit.load_model('mmt-files/kernik_2019_NaL_art.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    # Run for 20 s before running the VC protocol
    holding_proto = myokit.Protocol()
    for num in range(0, 50):
        holding_proto.add_step(hold_val, 100)
        holding_proto.add_step(hold_val + 1, 100)
    t_max = holding_proto.characteristic_time()
    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t_max, log_times=times)

    curr_names = ['ik1.i_K1',
                  'ito.i_to',
                  'ikr.i_Kr',
                  'iks.i_Ks',
                  'ical.i_CaL',
                  'icat.i_CaT',
                  'inak.i_NaK',
                  'ina.i_Na',
                  'inaca.i_NaCa',
                  'ipca.i_PCa',
                  'ifunny.i_f',
                  'ibna.i_b_Na',
                  'ibca.i_b_Ca',
                  'INaL.INaL',
                  'artifact.i_leak']


    tot_current = np.zeros(len(dat[curr_names[0]]))
    i_ion = np.zeros(len(dat[curr_names[0]]))

    for curr in curr_names:
        if curr == 'artifact.i_leak':
            continue
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array
            i_ion += np.array(dat[curr])

    if is_shown:
        idxs = [-1500, -500]
        diffs = [i_ion[idxs[1]] - i_ion[idxs[0]],
                 (dat['artifact.i_leak'][idxs[1]]
                                - dat['artifact.i_leak'][idxs[0]])/cm,
                 (dat['artifact.i_out'][idxs[1]]
                                - dat['artifact.i_out'][idxs[0]])/cm
                ]
        fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 8))
        fs = 16
        axs[0].plot(times, dat['artifact.v_cmd'])
        #axs[1].plot(times, tot_current, label="Abs(Total membrane currents)")
        axs[1].plot(times, i_ion, label=f'diff of {round(diffs[0], 2)}')
        axs[2].plot(times, np.array(dat['artifact.i_leak'])/cm,
                label=f'diff of {round(diffs[1], 2)}')
        axs[3].plot(times, np.array(dat['artifact.i_out'])/cm,
                label=f'diff of {round(diffs[2], 2)}')

        axs[0].set_ylabel('mV', fontsize=fs)
        axs[1].set_ylabel(r'$I_{ion}$ (A/F)', fontsize=fs)
        axs[2].set_ylabel(r'$I_{leak}$ (A/F)', fontsize=fs)
        axs[3].set_ylabel(r'$I_{ion}$ + $I_{leak}$ (A/F)', fontsize=fs)
        axs[3].set_xlabel('Time (ms)', fontsize=fs)

        for ax in axs:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.legend()

        #axs[1].legend()

        plt.show()




    curr_diff = np.abs(tot_current[98500] - tot_current[99500])

    if curr_diff < min_curr_diff:
        min_curr_diff = curr_diff
        best_holding = hold_val  
        #fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
        #fs = 16
        #axs[0].plot(times, dat['artifact.v_cmd'])
        #axs[0].plot(times, dat['membrane.V'])
        #axs[1].plot(times, tot_current, label="Abs(Total membrane currents)")
        #axs[1].plot(times,
        #        np.abs(np.array(dat['artifact.i_leak'])/cm), label="Abs(Seal leak)")

        #axs[0].set_ylabel('mV', fontsize=fs)
        #axs[1].set_ylabel('A/F', fontsize=fs)
        #axs[1].set_xlabel('Time (ms)', fontsize=fs)

        #for ax in axs:
        #    ax.spines['right'].set_visible(False)
        #    ax.spines['top'].set_visible(False)

        #axs[1].legend()

        #plt.show()

    artifact_curr_diff = np.abs(np.array(dat['artifact.i_leak'][98500]/cm) - 
                                    np.array(dat['artifact.i_leak'][99500]/cm))
    print(f'Holding value is {hold_val}, with average curr of {curr_diff}, and artifact current diff is {artifact_curr_diff}')


get_seal_comparison(-80)

#mod = myokit.load_model('mmt-files/kernik_2019_NaL_art.mmt')
#
#p = mod.get('engine.pace')
#p.set_binding(None)
#
#v_cmd = mod.get('artifact.v_cmd')
#v_cmd.set_rhs(0)
#v_cmd.set_binding('pace') # Bind to the pacing mechanism
#
## Run for 20 s before running the VC protocol
#holding_proto = myokit.Protocol()
#hold_val = -80
#holding_proto.add_step(-80, 30000)
#holding_proto.add_step(hold_val + 1, 500)
#t_max = holding_proto.characteristic_time()
#times = np.arange(0, t_max, 0.1)
#sim = myokit.Simulation(mod, holding_proto)
#dat = sim.run(t_max, log_times=times)
#
#curr_names = ['ik1.i_K1',
#              'ito.i_to',
#              'ikr.i_Kr',
#              'iks.i_Ks',
#              'ical.i_CaL',
#              'icat.i_CaT',
#              'inak.i_NaK',
#              'ina.i_Na',
#              'inaca.i_NaCa',
#              'ipca.i_PCa',
#              'ifunny.i_f',
#              'ibna.i_b_Na',
#              'ibca.i_b_Ca',
#              'INaL.INaL',
#              'artifact.i_leak']
#
#
#tot_current = np.zeros(len(dat[curr_names[0]]))
#
#for curr in curr_names:
#    if curr == 'artifact.i_leak':
#        continue
#    else:
#        as_array = np.abs(np.array(dat[curr]))
#        tot_current += as_array
#
#fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
#fs = 16
#        
#axs[0].plot(times, dat['artifact.v_cmd'])
#axs[1].plot(times, tot_current, label="Abs(Total membrane currents)")
#axs[1].plot(times,
#        np.abs(np.array(dat['artifact.i_leak'])/cm), label="Abs(Seal leak)")
#
#axs[0].set_ylabel('mV', fontsize=fs)
#axs[1].set_ylabel('A/F', fontsize=fs)
#axs[1].set_xlabel('Time (ms)', fontsize=fs)
#
#for ax in axs:
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)
#
#axs[1].legend()
#
#plt.show()

#
#
#
#fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
#piecewise_function, segment_dict, t_max = ind[0].get_myokit_protocol()
##t_max = 2000
#times = np.arange(0, t_max, 0.1)
#
#times = times[t_range[0]:t_range[1]]
#axs[0].plot(times, dat['artifact.v_cmd'][t_range[0]:t_range[1]], 'k--')
#axs[0].plot(times, dat['membrane.V'][t_range[0]:t_range[1]])
#axs[1].plot(times, np.array(
#    dat['artifact.i_out'][t_range[0]:t_range[1]]) / cm)
#
#currents = ['ina.i_Na', 'INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ito.i_to',
#                        'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks',
#                        'artifact.i_leak']
#
#for curr in curr_names:
#    if curr == 'artifact.i_leak':
#        contributions = np.abs(
#                np.array(dat[curr][t_range[0]:t_range[1]]) / cm) / (
#                        tot_current[t_range[0]:t_range[1]] )
#    else:
#        contributions = np.abs(
#                np.array(dat[curr][t_range[0]:t_range[1]])) / (
#                        tot_current[t_range[0]:t_range[1]])
#
#    if np.max(contributions) > .2:
#        axs[2].plot(times, contributions, label=curr, c=color_key[curr])
#    else:
#        axs[2].plot(times, contributions, 'k')
#
#if target_curr is not None:
#    times = np.arange(0, t_max, 0.1)
#    axs[2].axvspan(
#            times[max_idx-25], times[max_idx+25],facecolor='g', alpha=0.25)
#
#
#axs[0].set_ylabel('Vcmd')
#axs[1].set_ylabel('I_out')
#axs[2].set_ylabel('Contributions')
#axs[2].set_ylim(0, 1)
#axs[2].legend()
#
#if target_curr is not None:
#    fig.suptitle(
#            f'Max contribution for {target_curr} is {round(max_cont*100, 2)}%', fontsize=18)
#
#for ax in axs:
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)
#
#if is_shown:
#    plt.show()
#else:
#    return fig

