from utility import get_valid_cells
import pandas as pd
import numpy as np

'dmso_220829_001_3',# Great cell




def calc_summary_stats(drug_name):
    cells = get_valid_cells(drug_name)

    for f_name in cells:
        write_ap_meta(f_name)
        write_currs_dat(f_name)
        print(f_name)


### Utility functions
def write_ap_meta(f):
    ap_dat = pd.read_csv(f'./data/cells/{f}/ap_df.csv')
    ap_meta = pd.read_csv(f'./data/cells/{f}/ap_meta.csv')

    dvdt_is_na = np.array([get_dvdt(k, ap_dat[k]) for k in ap_meta['sweep'].values])
    all_mp = np.array([get_mdp(ap_dat[k]) for k in ap_meta['sweep'].values])
    all_apd90 = np.array([get_apd(ap_dat[k]) for k in ap_meta['sweep'].values])

    ap_meta['dvdt_max'] = dvdt_is_na[:, 0]
    ap_meta['is_upstroke'] = dvdt_is_na[:, 1]
    ap_meta['mdp'] = all_mp
    ap_meta['apd90'] = all_apd90[:, 0]

    ap_meta.to_csv(f'./data/cells/{f}/ap_summary_stats.csv', index=False)


def write_currs_dat(f):
    label_dict = {'Na1': 520, 'NaL': 2225, 'Na2': 2180, 'Kr': 3197, 'CaL': 3749, 'to': 4365, 'K1': 4955, 'f': 6388, 'Ks': 8905}

    curr_windows = {'Na1_min': [1000, 1050],
                    'Na2_min': [2660, 2700],
                    'NaL_avg': [2720, 2730],
                    'Kr_avg': [3690, 3697],
                    'CaL_min': [4248.5, 4260],
                    'to_max': [4862, 4882],
                    'K1_avg': [5452, 5460],
                    'f1_avg': [6870, 6885],
                    'f2_avg': [6889, 6895],
                    'Ks1_avg': [9350, 9403],
                    'Ks2_avg': [9405, 9410]}

    vc_dat = pd.read_csv(f'data/cells/{f}/vc_df.csv')
    vc_meta = pd.read_csv(f'data/cells/{f}/vc_meta.csv')

    curr_values = {}
    curr_values['sweep'] = vc_meta['sweep'].values
    curr_values['compound'] = vc_meta['compound'].values

    for curr_name, window in curr_windows.items():

        idxs = [int(val*25) for val in window]

        vc_slice = vc_dat.iloc[idxs[0]:idxs[1], :]
        if 'min' in curr_name:
            current_vals = vc_slice.min().values / vc_meta['cm'].values
        elif 'avg' in curr_name:
            current_vals = vc_slice.mean().values / vc_meta['cm'].values
        elif 'max' in curr_name:
            current_vals = vc_slice.max().values / vc_meta['cm'].values
        else:
            print('we have a problem')
            import pdb
            pdb.set_trace()

        curr_values[curr_name] = current_vals
        
    curr_values = pd.DataFrame(curr_values)
    curr_values.to_csv(f'./data/cells/{f}/vc_summary_stats.csv', index=False)

    return curr_values




def get_dvdt(ap_num, ap):
    stim_region = ap.values[int(48*25):int(65*25)]

    vals = moving_average(stim_region, 7)
    times = np.linspace(48,65, len(stim_region))
    times_avg = np.linspace(48, 65, len(vals))

    diff = np.diff(vals) / (times_avg[2]-times_avg[1])

    max_dvdt = np.max(diff)

    stim_start_idx = np.argmin(np.abs(times_avg - 50.5))
    stim_start_dvdt = np.mean(diff[stim_start_idx:(stim_start_idx+3)])

    max_dvdt -= stim_start_dvdt

    is_na_during_stim = not np.all(np.abs(np.diff(diff[7:20])) < 10)
    is_late_upstroke = not np.max(diff[21:35]) < 5

    if is_na_during_stim or is_late_upstroke:
        is_na = True
    else:
        is_na = False

    return [max_dvdt, is_na]


def get_apd(ap):
    stim_region = ap[int(48*25):int(65*25)]

    vals = moving_average(stim_region, 7)
    times = np.linspace(48,65, len(stim_region))
    times_avg = np.linspace(48, 65, len(vals))

    diff = np.diff(vals) / (times_avg[2]-times_avg[1])

    max_dvdt = np.max(diff)

    stim_start_idx = np.argmin(np.abs(times_avg - 50.5))
    stim_start_dvdt = np.mean(diff[stim_start_idx:(stim_start_idx+3)])

    max_dvdt -= stim_start_dvdt

    is_na_during_stim = not np.all(np.abs(np.diff(diff[7:20])) < 10)
    is_late_upstroke = not np.max(diff[21:35]) < 5

    if is_na_during_stim or is_late_upstroke:
        is_na = True
    else:
        is_na = False

    ap_orig = ap[:11000]
    times_orig = np.linspace(0, len(ap_orig)/25, len(ap_orig))

    ap = moving_average(ap[:11000], 4)

    max_v = np.max(ap)
    apa = max_v - np.min(ap)
    apd90_v = max_v - apa*.9
    argmax_v = np.argmax(ap)

    apd90_idx = np.argmin(np.abs(ap[argmax_v:] - apd90_v)) + argmax_v
    dvdt_max = np.argmax(np.diff(ap))

    #plt.plot(ap)
    #plt.axvline(dvdt_max, color='r')
    #plt.axvline(apd90_idx, color='b')

    times_new = np.linspace(0, 11000/25, len(ap))
    apd90 = times_new[apd90_idx] - times_new[dvdt_max]

    return [apd90, is_na]


def get_mdp(ap):
    return np.mean(ap[0:50])


def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]
    return np.array(new_vals)



def main():
    calc_summary_stats("dmso")
    calc_summary_stats("flecainide")
    calc_summary_stats("quinine")


if __name__ == "__main__":
    main()

