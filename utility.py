import myokit
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class VCProtocol():
    def __init__(self, segments):
        self.segments = segments #list of VCSegments

    def get_protocol_length(self):
        proto_length = 0
        for s in self.segments:
            proto_length += s.duration
        
        return proto_length

    def get_myokit_protocol(self):
        segment_dict = {'v0': '-82'}
        piecewise_txt = 'piecewise((engine.time >= 0 and engine.time < 500), v0, '
        current_time = 500

        #piecewise_txt = 'piecewise( '
        #current_time = 0
        #segment_dict = {}

        for i, segment in enumerate(self.segments):
            start = current_time
            end = current_time + segment.duration
            curr_step = f'v{i+1}'
            time_window = f'(engine.time >= {start} and engine.time < {end})'
            piecewise_txt += f'{time_window}, {curr_step}, '

            if segment.end_voltage is None:
                segment_dict[curr_step] = f'{segment.start_voltage}'
            else:
                slope = ((segment.end_voltage - segment.start_voltage) /
                                                                segment.duration)
                intercept = segment.start_voltage - slope * start

                segment_dict[curr_step] = f'{slope} * engine.time + {intercept}'
            
            current_time = end
        
        piecewise_txt += 'vp)'

        return piecewise_txt, segment_dict, current_time
        
    def plot_protocol(self, is_shown=False):
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        pts_v = []
        pts_t = []
        current_t = 0
        for seg in self.segments:
            pts_v.append(seg.start_voltage)
            if seg.end_voltage is None:
                pts_v.append(seg.start_voltage)
            else:
                pts_v.append(seg.end_voltage)
            pts_t.append(current_t)
            pts_t.append(current_t + seg.duration)

            current_t += seg.duration

        plt.plot(pts_t, pts_v)

        if is_shown:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xlabel('Time (ms)', fontsize=16)
            ax.set_xlabel('Voltage (mV)', fontsize=16)

            plt.show()

    def plot_with_curr(self, curr, cm=60):
        mod = myokit.load_model('mmt-files/kernik_2019_NaL_art.mmt')

        p = mod.get('engine.pace')
        p.set_binding(None)

        c_m = mod.get('artifact.c_m')
        c_m.set_rhs(cm)

        v_cmd = mod.get('artifact.v_cmd')
        v_cmd.set_rhs(0)
        v_cmd.set_binding('pace') # Bind to the pacing mechanism

        # Run for 20 s before running the VC protocol
        holding_proto = myokit.Protocol()
        holding_proto.add_step(-81, 30000)
        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)
        mod.set_state(sim.state())

        # Get protocol to run
        piecewise_function, segment_dict, t_max = self.get_myokit_protocol()
        mem = mod.get('artifact')

        for v_name, st in segment_dict.items():
            v_new = mem.add_variable(v_name)
            v_new.set_rhs(st)

        vp = mem.add_variable('vp')
        vp.set_rhs(0)

        v_cmd = mod.get('artifact.v_cmd')
        v_cmd.set_binding(None)
        vp.set_binding('pace')

        v_cmd.set_rhs(piecewise_function)
        times = np.arange(0, t_max, 0.1)
        ## CHANGE THIS FROM holding_proto TO SOMETHING ELSE
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t_max, log_times=times)

        fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
        axs[0].plot(times, dat['membrane.V'])
        axs[0].plot(times, dat['artifact.v_cmd'], 'k--')
        axs[1].plot(times, np.array(dat['artifact.i_out']) / cm)
        axs[2].plot(times, dat[curr])

        axs[0].set_ylabel('Voltage (mV)')
        axs[1].set_ylabel('I_out (A/F)')
        axs[2].set_ylabel(curr)
        for ax in axs:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
        plt.show()


class VCSegment():
    def __init__(self, duration, start_voltage, end_voltage=None):
        self.duration = duration
        self.start_voltage = start_voltage
        self.end_voltage = end_voltage


def get_valid_cells(drug_name):
    dmso_files = [
                  'dmso_220829_001_3',
                  'dmso_220829_003_1',
                  'dmso_220829_007_2',
                  'dmso_220929_001_1',
                  'dmso_230403_001_3',
                  'dmso_230403_003_1',
                  'dmso_230403_001_2'
                  ]

    flecainide_files = [
                    'flecainide_220617_001_3',
                    'flecainide_220617_005_2',
                    'flecainide_220830_010_3',
                    'flecainide_220901_009_3',
                    'flecainide_221028_003_4',
                    'flecainide_221028_008_4',
                    'flecainide_221117_001_2',
                    'flecainide_221117_001_3',
                    'flecainide_230403_005_4',
                    'flecainide_230403_009_3',
                    'flecainide_230403_007_3',
                    'flecainide_230403_008_1'
                    ]

    quinine_files = [
                    'quinine_220929_004_1',
                    'quinine_221013_008_2',
                    'quinine_221206_001_1',
                    'quinine_221206_001_3',
                    'quinine_221206_002_3',
                    'quinine_221206_003_1',
                    'quinine_221206_003_2',
                    'quinine_221206_007_2',
                    'quinine_221206_009_1',
                    'quinine_221206_009_2', 
                    'quinine_230120_007_3', 
                    'quinine_230120_010_2',
                    'quinine_230403_011_1',
                    'quinine_230403_011_3'
            ]

    if drug_name == 'dmso':
        return dmso_files
    if drug_name == 'flecainide':
        return flecainide_files
    if drug_name == 'quinine':
        return quinine_files


def get_vc_pts():
    steps = {
            'Voltage':  [-105, -36, -80, -112, -2, -116, -31, -80, 24, -37, 28, -80, 3, -90, -80, -68, 57, -80, -68, -120, -80, -120, -120, -77, -80, 23, 40, -1, -80, -110, -20, -80, -110, .1, -80, -110, 20, -80, -45, -10, -80, -45, 5, -80, -45, 20, -80, 30, -50, 10, -80, 30, -35, 10, -80, 30, -20, 10, -80],
            'Duration': [12, 40, 500, 5, 364, 754, 88, 500, 422, 14, 40, 500, 10, 40, 500, 71, 43, 500, 45, 46, 500, 18, 876, 40, 500, 357, 1620, 40, 500, 30, 50, 500, 30, 50, 500, 30, 50, 500, 400, 100, 500, 400, 100, 500, 400, 100, 500, 300, 10, 50, 500, 300, 10, 50, 500, 300, 10, 50, 20]
            }

    #steps['Voltage'] = steps['Voltage'][0:29]
    #steps['Duration'] = steps['Duration'][0:29]

    print(pd.DataFrame(steps))
    
    #import pdb
    #pdb.set_trace()

    pts = [[0, -80]]
    t = 999 
    curr_voltage = -80

    for i in range(0, len(steps['Voltage'])):
        pts.append([t, curr_voltage])
        curr_voltage = int(float(steps['Voltage'][i]))
        pts.append([t, curr_voltage])

        t += steps['Duration'][i] 

    return np.array(pts)


def get_model_response(model_name, scale_params=None):
    #paci_names = 
    #import pdb
    #pdb.set_trace()
    if model_name == 'Kernik':
        mod = myokit.load_model(
                './vc_proto_optimization/mmt-files/kernik_2019_NaL_art.mmt')
        mod['geom']['Cm'].set_rhs(20)
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
    else:
        mod = myokit.load_model(
                './vc_proto_optimization/mmt-files/paci_cardio_lei.mmt')

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

    return dat
 
