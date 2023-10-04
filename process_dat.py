import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import heka_reader
from os import mkdir, listdir


def save_cell_data(folder, date, file_num, ch, drug):
    # vc_plot.pdf
    # ap_plot.pdf
    # vc_dat.csv
    # ap_dat.csv
    # meta.csv
    f_path = f'{folder}/{date}_{file_num}.dat'
    ap_meta_path = f'{folder}/{date}/{date}_{file_num}_1CCstim.xls'
    vc_meta_path = f'{folder}/{date}/{date}_{file_num}_1optTest.xls'
    ap_meta = pd.read_csv(ap_meta_path, sep=r'\t')
    #try:
    with open(vc_meta_path) as f:
        lines = f.read()
        first = lines.split('\n', 1)[0]

    names = first.split('\t')

    if len(names) < 38:
        names += [f'Online6{n}' for n in [1, 2, 3, 4]]

    vc_meta = pd.read_csv(vc_meta_path, sep=r'\t', names=names)

    if ap_meta[f'Age({ch})'].isnull().values.any():
        print('CHANNEL IS NULL')

    bundle = heka_reader.Bundle(f_path)
    
    ap_df, ap_meta = get_ap_data(bundle, ap_meta, ch)
    vc_df, vc_meta = get_vc_data(bundle, vc_meta, ch)

    all_csv_cells = listdir(f'data/cells')
    new_folder = f'data/cells/{drug}_{date}_{file_num}_{ch}'

    if new_folder.split('/')[-1] not in all_csv_cells:
        mkdir(new_folder)

    ap_df.to_csv(f'{new_folder}/ap_df.csv', index=False)
    ap_meta.to_csv(f'{new_folder}/ap_meta.csv', index=False)
    vc_df.to_csv(f'{new_folder}/vc_df.csv', index=False)
    vc_meta.to_csv(f'{new_folder}/vc_meta.csv', index=False)


def get_ap_data(bundle, ap_meta, ch):
    ap_indices = {}

    #iterate through all trials
    labels = {}

    for i, trial in enumerate(bundle.pul[0]):
        labels[f'{trial.Label}_{i}'] = [trial.SeriesCount, trial.NumberSweeps]

    keys = [k for k in labels.keys() if 'Rack' in k]

    if keys:
        all_keys = []
        for i in range(1, 5):
            curr_keys = [k for k in keys if f'C:{i}' in k]

            if curr_keys:
                all_keys.append(curr_keys[-1])
    else:
        all_keys = [k for k in labels.keys() if 'M' in k]

    all_meta = {'sweep': [], 'time': [], 'compound': [], 'age': []}
    all_dat = {}

    for curr_val, k in enumerate(all_keys):
        v = labels[k]
        rows = ap_meta[ap_meta['Sweep'].str.contains(f'1_{v[0]}_')]

        sweeps = rows['Sweep'].values

        for sw_range in [sweeps[0:12], sweeps[-12:]]:
            for sweep in sw_range:
                ins = [int(val) for val in sweep.split('_')]
                dat = bundle.data[ins[0]-1, ins[1]-1, ins[2]-2, ch-1]
                dat = dat * 1000

                #all_meta: Compound, Online1 (), Online3 (), Online 4 (), Online 7 (GK1 I think), Online8 ()
                sweep_meta = ap_meta[ap_meta['Sweep'] == sweep]
                all_meta['sweep'].append(sweep_meta['Sweep'].values[0])
                all_meta['time'].append(sweep_meta['Time'].values[0])
                if (('wash' in sweep_meta[f'Compound({ch})'].values[0]) or
                        ('F:3' in sweep_meta[f'Compound({ch})'].values[0])):
                    all_meta['compound'].append('Baseline')
                elif 'C:1' in sweep_meta[f'Compound({ch})'].values[0]:
                    all_meta['compound'].append('3xEFPC')
                elif 'C:2' in sweep_meta[f'Compound({ch})'].values[0]:
                    all_meta['compound'].append('10xEFPC')
                elif 'C:3' in sweep_meta[f'Compound({ch})'].values[0]:
                    all_meta['compound'].append('20xEFPC')
                elif 'C:4' in sweep_meta[f'Compound({ch})'].values[0]:
                    all_meta['compound'].append('Washout')
                else:
                    print('MISSED SOMETHING')
                    import pdb
                    pdb.set_trace()

                all_meta['age'].append(sweep_meta[f'Age({ch})'].values[0])

                all_dat[sweep_meta['Sweep'].values[0]] = dat

    dat_df = pd.DataFrame(all_dat)
    meta_df = pd.DataFrame(all_meta)

    return dat_df, meta_df


def get_vc_data(bundle, vc_meta, ch):
    vc_indices = {}

    #iterate through all trials
    labels = {}

    for i, trial in enumerate(bundle.pul[0]):
        labels[f'{trial.Label}_{i}'] = [trial.SeriesCount, trial.NumberSweeps]

    keys = [k for k in labels.keys() if 'optIPSC' in k]

    all_meta = {'sweep': [], 'time': [], 'compound': [], 'age': [], 'cm': [], 'r_series': [], 'r_seal': []}
    all_dat = {}

    for index, k in enumerate(keys):
        v = labels[k]
        rows = vc_meta[vc_meta['Sweep'].str.contains(f'1_{v[0]}_')]

        sweeps = rows['Sweep'].values

        for sweep in sweeps:
            ins = [int(val) for val in sweep.split('_')]
            dat = bundle.data[ins[0]-1, ins[1]-1, ins[2]-2, ch-1]
            dat = dat

            sweep_meta = vc_meta[vc_meta['Sweep'] == sweep]
            all_meta['sweep'].append(sweep_meta['Sweep'].values[0])
            all_meta['time'].append(sweep_meta['Time'].values[0])

            if (('wash' in sweep_meta[f'Compound({ch})'].values[0]) or
                        ('F:3' in sweep_meta[f'Compound({ch})'].values[0])):
                all_meta['compound'].append('Baseline')
            elif 'C:1' in sweep_meta[f'Compound({ch})'].values[0]:
                all_meta['compound'].append('3xEFPC')
            elif 'C:2' in sweep_meta[f'Compound({ch})'].values[0]:
                all_meta['compound'].append('10xEFPC')
            elif 'C:3' in sweep_meta[f'Compound({ch})'].values[0]:
                all_meta['compound'].append('20xEFPC')
            else:
                all_meta['compound'].append('Washout')

            all_meta['age'].append(sweep_meta[f'Age({ch})'].values[0])
            all_meta['cm'].append(sweep_meta[f'Online1({ch})'].values[0])
            all_meta['r_series'].append(sweep_meta[f'Online2({ch})'].values[0])
            all_meta['r_seal'].append(sweep_meta[f'Online4({ch})'].values[0])

            all_dat[sweep_meta['Sweep'].values[0]] = dat

    dat_df = pd.DataFrame(all_dat)
    meta_df = pd.DataFrame(all_meta)

    return dat_df, meta_df


def write_dmso():
    #8/29
    #save_cell_data(folder='./data/iPSC-dmso_8-29-2022', date='220829', file_num='001', ch=1, drug='dmso')
    #save_cell_data(folder='./data/iPSC-dmso_8-29-2022', date='220829', file_num='001', ch=2, drug='dmso')
    #save_cell_data(folder='./data/iPSC-dmso_8-29-2022', date='220829', file_num='001', ch=3, drug='dmso')

    #save_cell_data(folder='./data/iPSC-dmso_8-29-2022', date='220829', file_num='002', ch=1, drug='dmso')

    #save_cell_data(folder='./data/iPSC-dmso_8-29-2022', date='220829', file_num='003', ch=1, drug='dmso')
    #save_cell_data(folder='./data/iPSC-dmso_8-29-2022', date='220829', file_num='003', ch=3, drug='dmso')

    #save_cell_data(folder='./data/iPSC-dmso_8-29-2022', date='220829', file_num='007', ch=2, drug='dmso')

    ##9/29 
    #save_cell_data(folder='./data/iPSC-quinine_9-29-2022', date='220929', file_num='001', ch=1, drug='dmso')
    #save_cell_data(folder='./data/iPSC-quinine_9-29-2022', date='220929', file_num='001', ch=2, drug='dmso')
    #save_cell_data(folder='./data/iPSC-quinine_9-29-2022', date='220929', file_num='001', ch=3, drug='dmso')

    #save_cell_data(folder='./data/iPSC-quinine_9-29-2022', date='220929', file_num='002', ch=2, drug='dmso')
    #save_cell_data(folder='./data/iPSC-quinine_9-29-2022', date='220929', file_num='002', ch=3, drug='dmso')

    #save_cell_data(folder='./data/iPSC-quinine_9-29-2022', date='220929', file_num='003', ch=1, drug='dmso')

    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='002', ch=4, drug='dmso')

    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='003', ch=2, drug='dmso')
    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='003', ch=4, drug='dmso')

    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='004', ch=3, drug='dmso')
    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='004', ch=4, drug='dmso')
    #
    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='007', ch=3, drug='dmso')
    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='007', ch=4, drug='dmso')

    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='009', ch=2, drug='dmso')
    #save_cell_data(folder='./data/iPSC-dmso_11-15-2022', date='221115', file_num='009', ch=3, drug='dmso')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='001', ch=1, drug='dmso')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='001', ch=2, drug='dmso')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='001', ch=3, drug='dmso')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='003', ch=1, drug='dmso')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='003', ch=3, drug='dmso')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='003', ch=4, drug='dmso')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='004', ch=1, drug='dmso')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='004', ch=3, drug='dmso')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='004', ch=4, drug='dmso')

    save_cell_data(folder='./data/iPSC_dmso_4-12-2023', date='230412', file_num='005', ch=2, drug='dmso')




def write_quinine():
    #9/29
    #save_cell_data(folder='./data/IPSC-quinine_9-29-2022', date='220929', file_num='004', ch=1, drug='quinine')
    #save_cell_data(folder='./data/IPSC-quinine_9-29-2022', date='220929', file_num='004', ch=2, drug='quinine')

    #save_cell_data(folder='./data/IPSC-quinine_9-29-2022', date='220929', file_num='005', ch=2, drug='quinine')
    #save_cell_data(folder='./data/IPSC-quinine_9-29-2022', date='220929', file_num='005', ch=3, drug='quinine')

    ##10/13
    #save_cell_data(folder='./data/iPSC-quinine_10-13-2022', date='221013', file_num='002', ch=1, drug='quinine')
    #save_cell_data(folder='./data/iPSC-quinine_10-13-2022', date='221013', file_num='002', ch=3, drug='quinine')

    #save_cell_data(folder='./data/iPSC-quinine_10-13-2022', date='221013', file_num='006', ch=2, drug='quinine')
    #
    #save_cell_data(folder='./data/iPSC-quinine_10-13-2022', date='221013', file_num='007', ch=2, drug='quinine')

    #save_cell_data(folder='./data/iPSC-quinine_10-13-2022', date='221013', file_num='008', ch=1, drug='quinine')
    #save_cell_data(folder='./data/iPSC-quinine_10-13-2022', date='221013', file_num='008', ch=2, drug='quinine')

    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='001', ch=1, drug='quinine')
    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='001', ch=2, drug='quinine')
    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='001', ch=3, drug='quinine')
    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='001', ch=4, drug='quinine')

    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='002', ch=3, drug='quinine')

    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='003', ch=1, drug='quinine')
    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='003', ch=2, drug='quinine')

    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='007', ch=2, drug='quinine')

    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='009', ch=1, drug='quinine')
    #save_cell_data(folder='./data/iPSC-quinine_12-6-2022', date='221206', file_num='009', ch=2, drug='quinine')


   # save_cell_data(folder='./data/iPSC-quinine_1-20-2023', date='230120', file_num='003', ch=3, drug='quinine')

   # save_cell_data(folder='./data/iPSC-quinine_1-20-2023', date='230120', file_num='006', ch=1, drug='quinine')
   # save_cell_data(folder='./data/iPSC-quinine_1-20-2023', date='230120', file_num='006', ch=3, drug='quinine')

   # save_cell_data(folder='./data/iPSC-quinine_1-20-2023', date='230120', file_num='007', ch=3, drug='quinine')

   # save_cell_data(folder='./data/iPSC-quinine_1-20-2023', date='230120', file_num='008', ch=4, drug='quinine')

   # save_cell_data(folder='./data/iPSC-quinine_1-20-2023', date='230120', file_num='010', ch=1, drug='quinine')
   # save_cell_data(folder='./data/iPSC-quinine_1-20-2023', date='230120', file_num='010', ch=2, drug='quinine')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='011', ch=1, drug='quinine')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='011', ch=3, drug='quinine')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='012', ch=2, drug='quinine')


def write_flecainide():
    #6/17
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='001', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='001', ch=3, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='001', ch=4, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='002', ch=2, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='003', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='003', ch=2, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='003', ch=4, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='004', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='004', ch=2, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='005', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='005', ch=2, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='005', ch=3, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='006', ch=2, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='006', ch=3, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_6-17-2022', date='220617', file_num='006', ch=4, drug='flecainide')

    ##8/30
    #save_cell_data(folder='./data/iPSC-flecainide_8-30-2022', date='220830', file_num='007', ch=3, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_8-30-2022', date='220830', file_num='010', ch=3, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_8-30-2022', date='220830', file_num='012', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_8-30-2022', date='220830', file_num='012', ch=2, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_8-30-2022', date='220830', file_num='012', ch=3, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_8-30-2022', date='220830', file_num='012', ch=4, drug='flecainide')

    ##9/1
    #save_cell_data(folder='./data/iPSC-flecainide_9-01-2022', date='220901', file_num='009', ch=3, drug='flecainide')

    ##10/29
    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='003', ch=3, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='003', ch=4, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='006', ch=2, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='006', ch=3, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='007', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='007', ch=3, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='007', ch=4, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='008', ch=3, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_10-28-2022', date='221028', file_num='008', ch=4, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='001', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='001', ch=2, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='001', ch=3, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='002', ch=1, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='003', ch=1, drug='flecainide')
    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='003', ch=2, drug='flecainide')
    ##save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='003', ch=3, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='005', ch=4, drug='flecainide')

    #save_cell_data(folder='./data/iPSC-flecainide_11-17-2022', date='221117', file_num='007', ch=2, drug='flecainide')


    #HCM from 2/15/2022
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-15-2023', date='230215', file_num='001', ch=1, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-15-2023', date='230215', file_num='001', ch=2, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-15-2023', date='230215', file_num='001', ch=3, drug='hcm_flecainide')

    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-15-2023', date='230215', file_num='002', ch=2, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-15-2023', date='230215', file_num='002', ch=3, drug='hcm_flecainide')

    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-15-2023', date='230215', file_num='003', ch=1, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-15-2023', date='230215', file_num='003', ch=3, drug='hcm_flecainide')


    #HCM from 2/17/2022
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='001', ch=1, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='001', ch=2, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='001', ch=3, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='001', ch=4, drug='hcm_flecainide')

    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='003', ch=1, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='003', ch=3, drug='hcm_flecainide')

    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='004', ch=1, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='004', ch=2, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='004', ch=4, drug='hcm_flecainide')

    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='005', ch=1, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='005', ch=2, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='005', ch=3, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='005', ch=4, drug='hcm_flecainide')

    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='006', ch=1, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='006', ch=3, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='006', ch=4, drug='hcm_flecainide')

    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='007', ch=3, drug='hcm_flecainide')
    #save_cell_data(folder='./data/iPSC_HCM_flecainide_2-17-2023', date='230217', file_num='007', ch=4, drug='hcm_flecainide')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='005', ch=1, drug='flecainide')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='005', ch=3, drug='flecainide')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='005', ch=4, drug='flecainide')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='006', ch=2, drug='flecainide')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='006', ch=3, drug='flecainide')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='007', ch=3, drug='flecainide')
    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='007', ch=4, drug='flecainide')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='008', ch=1, drug='flecainide')

    save_cell_data(folder='./data/iPSC-dmso-flec-quin_4-3-2023', date='230403', file_num='009', ch=3, drug='flecainide')


def write_moxifloxacin():
    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='001', ch=2, drug='moxifloxacin')
    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='001', ch=4, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='003', ch=2, drug='moxifloxacin')
    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='003', ch=3, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='004', ch=1, drug='moxifloxacin')
    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='004', ch=2, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='006', ch=2, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='007', ch=1, drug='moxifloxacin')
    save_cell_data(folder='./data/iPSC-moxifloxacin_6-22-2022', date='220622', file_num='007', ch=4, drug='moxifloxacin')

    #8/31
    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='001', ch=3, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='003', ch=1, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='006', ch=1, drug='moxifloxacin')
    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='006', ch=3, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='008', ch=1, drug='moxifloxacin')
    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='008', ch=3, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='011', ch=1, drug='moxifloxacin')

    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='014', ch=1, drug='moxifloxacin')
    save_cell_data(folder='./data/iPSC-moxifloxicin_8-31-2022', date='220831', file_num='014', ch=4, drug='moxifloxacin')


def write_verapamil():
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='002', ch=4, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='003', ch=2, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='003', ch=3, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='005', ch=1, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='005', ch=3, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='006', ch=3, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='007', ch=1, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='007', ch=2, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='007', ch=3, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='007', ch=4, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='008', ch=1, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='008', ch=4, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='010', ch=1, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='011', ch=1, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='011', ch=2, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='011', ch=3, drug='verapamil')

    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='012', ch=3, drug='verapamil')
    save_cell_data(folder='./data/iPSC-verapamil_6-21-2022', date='220621', file_num='012', ch=4, drug='verapamil')


def write_chloroquine():
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='001', ch=1, drug='chloroquine')
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='001', ch=2, drug='chloroquine')
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='001', ch=3, drug='chloroquine')
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='001', ch=4, drug='chloroquine')

    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='002', ch=2, drug='chloroquine')

    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='005', ch=1, drug='chloroquine')
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='005', ch=2, drug='chloroquine')
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='005', ch=4, drug='chloroquine')

    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='006', ch=3, drug='chloroquine')

    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='007', ch=1, drug='chloroquine')
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='007', ch=3, drug='chloroquine')

    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='009', ch=4, drug='chloroquine')

    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='010', ch=2, drug='chloroquine')
    save_cell_data(folder='./data/iPSC-chloroquine-6-24-2022', date='220624', file_num='010', ch=4, drug='chloroquine')


#write_dmso()
#write_moxifloxacin()
write_quinine()
#write_flecainide()
#write_quinine()
#write_chloroquine()
#write_verapamil()
