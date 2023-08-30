import pickle
from deap import base, creator, tools
from vc_opt_ga import simulate_model
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from utility_classes import VCProtocol, VCSegment
from os import listdir



class MyCustomUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try:
            creator.create('FitnessMax', base.Fitness, weights=(1.0,))
            creator.create('Individual', list, fitness=creator.FitnessMax)
            if module == '__main__':
                module = 'vc_opt_ga'
            return super().find_class(module, name)
        except:
            import pdb
            pdb.set_trace()


def plot_with_curr_contribution(ind, target_curr, is_shown=False, cm=60):
    cm=60
    dat = simulate_model(ind)

    if dat is None:
        print("Error in model run")
        return 0

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
    for curr in curr_names:
        if curr == 'artifact.i_leak':
            as_array = np.abs(np.array(dat[curr])) / cm
            tot_current += as_array
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array


    contributions = np.abs(np.array(dat[target_curr])) / tot_current

    max_contrib = np.max(contributions)
    max_arg = np.argmax(contributions)

    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 8))
    fig.suptitle(f'{round(max_contrib*100, 2)}% at {max_arg/10} ms', fontsize=16)
    piecewise_function, segment_dict, t_max = ind[0].get_myokit_protocol()
    #t_max = 2000
    times = np.arange(0, t_max, 0.1)

    axs[0].plot(times, dat['artifact.v_cmd'], 'k--')
    axs[0].plot(times, dat['membrane.V'])
    axs[1].plot(times, np.array(dat['artifact.i_out']) / cm)
    axs[2].plot(times, np.array(dat[target_curr]))
    axs[3].plot(times, contributions)

    axs[0].set_ylabel('Vcmd')
    axs[1].set_ylabel('I_out')
    axs[2].set_ylabel(target_curr)
    axs[3].set_ylabel('Contribution')
    axs[3].set_ylim(0, 1)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    if is_shown:
        plt.show()
    else:
        return fig


def get_best_proto(folder, f_name):
    with open(f'results/{folder}/{f_name}', 'rb') as f:
        unpickler = MyCustomUnpickler(f)
        obj = unpickler.load()

    pop = pickle.load(open(f'results/{folder}/{f_name}', "rb" ))
    last_gen = pop[-1]
    last_gen.sort(key=lambda x: x.fitness.values[0], reverse=True)

    return last_gen[0]


def plot_vc_protocols(folder='exp_2'):
    files = ['INaL.pkl', 'ikr.pkl', 'ical.pkl', 'ina.pkl', 'ito.pkl', 'ik1.pkl', 'ifunny.pkl', 'iks.pkl']
    currents = ['INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ina.i_Na', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks'] 

    for f_name, target_curr in dict(zip(files, currents)).items():
        if f_name not in listdir(f'results/{folder}'):
            continue
        ind = get_best_proto(folder, f_name)

        plot_with_curr_contribution(ind, target_curr=target_curr)

        plt.savefig(f'results/{folder}/{f_name.split(".")[0]}.pdf')


def shorten_vc_protocols(folder):
    files = ['ina.pkl', 'INaL.pkl', 'ikr.pkl', 'ical.pkl',                             'ito.pkl', 'ik1.pkl', 'ifunny.pkl', 'iks.pkl']
    currents = ['ina.i_Na', 'INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks'] 

    new_protos = {}

    for f_name, target_curr in dict(zip(files, currents)).items():
        ind = get_best_proto(folder, f_name)
        max_iso, max_idx = get_current_contribution(ind[0], target_curr)
        print(f'Start max for {target_curr}: {max_iso}')
        new_proto = remove_end_of_proto(ind[0], target_curr)
        shortened_proto = remove_start_of_proto(new_proto, target_curr)

        max_iso, max_idx = get_current_contribution(shortened_proto,
                                                        target_curr)
        print(f'End max for {target_curr}: {max_iso}')

        new_protos[target_curr] = shortened_proto
        #shortened_proto.plot_with_curr(target_curr)
        #plot_with_curr_contribution([shortened_proto],target_curr,is_shown=True)


    pickle.dump(new_protos, open(f"results/{folder}/shortened_protocols.pkl", "wb" ))

    return new_protos 


def remove_end_of_proto(proto, target_curr):
    max_cont, max_idx = get_current_contribution(proto, target_curr)
    # The max_idx was calculated with a 500 ms -80 mV holding step before it
    max_idx -= 5000

    # Remove everything in protocol from end of max_idx to end
    # Systematically shorten the front of the protocol
    tot_time, last_seg_start = 0, 0
    max_proto_time = (max_idx / 10)
    end_proto_time = max_proto_time + 40
    new_proto = VCProtocol(segments=[])
    for which_seg, seg in enumerate(proto.segments):
        tot_time += seg.duration
        if tot_time > end_proto_time:
            dur = end_proto_time - last_seg_start
            new_proto.segments.append(VCSegment(dur,
                                        seg.start_voltage, seg.end_voltage))
            break
        else:
            new_proto.segments.append(seg)

        last_seg_start = tot_time

    return new_proto


def remove_start_of_proto(proto, target_curr, reduction_threshold=.05):
    trim_size = 50
    continue_shortening = True
    max_isolation, max_idx = get_current_contribution(proto, target_curr)


    print(f'Protocol length is: {proto.get_protocol_length()}')
    #print(f'\tProtocol contribution is: {max_isolation}')

    while continue_shortening is True:
        new_proto = trim_proto_start(proto, trim_size=trim_size) 

        new_proto_isolation, new_idx = get_current_contribution(new_proto,
                                                                    target_curr)

        max_time = new_idx / 10 - 500
        if max_time < trim_size:
            return proto


        if (max_isolation - new_proto_isolation) > (
                                    max_isolation * reduction_threshold):
            continue_shortening = False
            continue

        if new_proto.get_protocol_length() < trim_size:
            continue_shortening = False
            continue

        proto = new_proto
        print(f'New protocol length is: {proto.get_protocol_length()}')
        #print(f'\tNew protocol contribution is: {new_proto_isolation}')

    return proto


def trim_proto_start(proto, trim_size=50):
    seg_end_points = [s.duration for s in proto.segments]
    accum = np.array([sum(seg_end_points[0:i])
                for i, j in enumerate(proto.segments)]) + seg_end_points[0]
    which_seg = np.argmax((accum-trim_size)>0)

    segments = []
    for s in proto.segments[which_seg:]:
        segments.append(VCSegment(s.duration, s.start_voltage, s.end_voltage))

    segments[0].duration = accum[which_seg] - trim_size
    new_proto = VCProtocol(segments)

    return new_proto


def create_long_proto(all_protos, folder):
    all_segs = []
    for curr, proto in all_protos.items():
        for s in proto.segments:
            all_segs.append(VCSegment(
                            s.duration, s.start_voltage, s.end_voltage))
        all_segs.append(VCSegment(500, -80))

    long_proto = VCProtocol(all_segs)
    long_proto.plot_protocol(True)
    pickle.dump(long_proto, open(f"results/{folder}/long_proto.pkl", "wb" ))


def compare_contributions(folder):
    files = ['ina.pkl', 'INaL.pkl', 'ikr.pkl', 'ical.pkl',                             'ito.pkl', 'ik1.pkl', 'ifunny.pkl', 'iks.pkl']
    currents = ['ina.i_Na', 'INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks'] 

    long_proto = pickle.load(open(f'./results/{folder}/long_proto.pkl', 'rb'))

    for f_name, target_curr in dict(zip(files, currents)).items():
        ind = get_best_proto(folder, f_name)


        max_isolation, max_idx = get_current_contribution(ind[0], target_curr)
        long_max, long_idx = get_current_contribution(long_proto, target_curr)

        print(
          f'{target_curr}: Original max: {max_isolation} – long max: {long_max} at {long_idx/10}')
    

def plot_cal_curr(folder):
    current = 'ical.i_CaL' 

    long_proto = pickle.load(open(f'./results/{folder}/long_proto.pkl', 'rb'))

    dat = simulate_model([long_proto])

    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['artifact.v_cmd'])
    axs[1].plot(dat['ical.i_CaL'])
    axs[2].plot(dat['ical.d'])
    axs[3].plot(dat['ical.f'])

    axs[0].set_ylabel('mV')
    axs[1].set_ylabel('Current')
    axs[2].set_ylabel('d')
    axs[3].set_ylabel('f')
    plt.show()


def plot_proto(ind, is_shown=True):
    cm=60
    dat = simulate_model(ind)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    piecewise_function, segment_dict, t_max = ind[0].get_myokit_protocol()
    times = np.arange(0, t_max, 0.1)

    axs[0].plot(times, dat['artifact.v_cmd'], 'k--')
    axs[0].plot(times, dat['membrane.V'])
    axs[1].plot(times, np.array(dat['artifact.i_out']) / cm)

    fs = 16
    axs[0].set_ylabel('mV', fontsize=fs)
    axs[1].set_ylabel('I_out (A/F)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    if is_shown:
        plt.show()
    else:
        return fig


def get_current_contribution(proto, target_current, window=2):
    dat = simulate_model([proto])
    capacitance = 60

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
    for curr in curr_names:
        if curr == 'artifact.i_leak':
            as_array = np.abs(np.array(dat[curr])) / capacitance
            tot_current += as_array
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array


    contrib = np.abs(np.array(dat[target_current])) / tot_current

    if window > 1:
        smoothed_contrib = (np.convolve(contrib,
                                np.ones(window), 'valid') / window)
        max_isolation = np.max(smoothed_contrib)
        max_idx = np.argmax(smoothed_contrib)

    else:
        max_isolation = np.max(contrib)
        max_idx = np.argmax(contrib)

    return max_isolation, max_idx 


def plot_all_curr_contributions(ind, target_curr=None, is_shown=True):
    cm=60
    dat = simulate_model(ind)

    if dat is None:
        print("Error in model run")
        return 0
    
    color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'brown', 'lightsalmon',
                    'goldenrod', 'greenyellow', 'aquamarine', 'turquoise',
                    'lightsteelblue', 'rebeccapurple', 'fuchsia'] 

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

    color_key = dict(zip(curr_names, color_list))

    tot_current = np.zeros(len(dat[curr_names[0]]))

    if target_curr is None:
        t_range = [0, len(dat[curr_names[0]])]
    else:
        max_cont, max_idx= get_current_contribution(ind[0], target_curr)

        t_range = [max_idx- 1000, max_idx + 1000]
    
    for curr in curr_names:
        if curr == 'artifact.i_leak':
            as_array = np.abs(np.array(dat[curr])) / cm
            tot_current += as_array
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array


    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    piecewise_function, segment_dict, t_max = ind[0].get_myokit_protocol()
    #t_max = 2000
    times = np.arange(0, t_max, 0.1)

    times = times[t_range[0]:t_range[1]]
    axs[0].plot(times, dat['artifact.v_cmd'][t_range[0]:t_range[1]], 'k--')
    axs[0].plot(times, dat['membrane.V'][t_range[0]:t_range[1]])
    axs[1].plot(times, np.array(
        dat['artifact.i_out'][t_range[0]:t_range[1]]) / cm)

    currents = ['ina.i_Na', 'INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks',
                            'artifact.i_leak'] 

    for curr in curr_names:
        if curr == 'artifact.i_leak':
            contributions = np.abs(
                    np.array(dat[curr][t_range[0]:t_range[1]]) / cm) / (
                            tot_current[t_range[0]:t_range[1]] )
        else:
            contributions = np.abs(
                    np.array(dat[curr][t_range[0]:t_range[1]])) / (
                            tot_current[t_range[0]:t_range[1]])

        if np.max(contributions) > .2:
            axs[2].plot(times, contributions, label=curr, c=color_key[curr])
        else:
            axs[2].plot(times, contributions, 'k')

    if target_curr is not None:
        times = np.arange(0, t_max, 0.1)
        axs[2].axvspan(
                times[max_idx-25], times[max_idx+25],facecolor='g', alpha=0.25)

    axs[0].set_ylabel('Vcmd')
    axs[1].set_ylabel('I_out')
    axs[2].set_ylabel('Contributions')
    axs[2].set_ylim(0, 1)
    axs[2].legend()

    if target_curr is not None:
        fig.suptitle(
                f'Max contribution for {target_curr} is {round(max_cont*100, 2)}%', fontsize=18)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    if is_shown:
        plt.show()
    else:
        return fig


def make_contributions_pdf(folder):
    path = f'./results/{folder}'
    long_proto = pickle.load(open(f'{path}/long_proto.pkl', 'rb'))
    currents = ['ina.i_Na', 'INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks'] 

    for seg in long_proto.segments:
        print(f'Duration: {round(seg.duration, 1)}, Voltage: {round(seg.start_voltage, 1)}')

    all_figs = []
    all_figs.append(plot_proto([long_proto], is_shown=False))



    for curr in currents:
        fig = plot_all_curr_contributions([long_proto], curr, is_shown=False)
        all_figs.append(fig)

    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages(
                                        f"{path}/contributions_{folder}.pdf")

    for fig in all_figs:
        pdf.savefig(fig)
    pdf.close()


def proto_steps_to_csv(folder):
    long_proto = pickle.load(open(f'results/{folder}/long_proto.pkl', 'rb'))

    f = open(f'results/{folder}/long_proto_steps.csv', "a")

    for i, seg in enumerate(long_proto.segments):
        f.write(f'Step {i+2}, Voltage: {round(seg.start_voltage, 1)}, Duration: {round(seg.duration, 1)}\n')

    f.close()

    
folder = 'exp_14' 
#all_protos = shorten_vc_protocols(folder=folder)
#create_long_proto(all_protos, folder=folder)
long_proto = pickle.load(open(f'results/{folder}/long_proto.pkl', 'rb'))
#proto_steps_to_csv(folder)
#plot_vc_protocols(folder)
#compare_contributions(folder)
#plot_proto([long_proto])
#plot_all_curr_contributions([long_proto], is_shown=True)
make_contributions_pdf(folder)

#plot_cal_curr(folder)
