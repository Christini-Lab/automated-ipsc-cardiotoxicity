import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from seaborn import pointplot, swarmplot 
from scipy.stats import ttest_ind


from utility import get_valid_cells
import matplotlib

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10 
plt.rcParams['axes.labelsize'] = 10 
plt.rc('legend', fontsize = 8)
matplotlib.rcParams['pdf.fonttype'] = 42



def plot_figure(feature, curr_name):
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6))
    fig.subplots_adjust(.1, .15, .95, .95)

    axs = axs.flatten()


    plot_ap_feature(axs[0], 'flecainide', feature, is_pct_change=True)
    plot_ap_feature(axs[1], 'quinine', feature, is_pct_change=True)

    window = [2674.3, 2685]
    channel = 'I_Na'

    plot_curr(axs[2], 'flecainide', window, channel, is_change=True)
    plot_curr(axs[3], 'quinine', window, channel, is_change=True)

    for vals in [[0, 'A'], [1, 'B'], [2, 'C'], [3, 'D']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)
    
    axs[0].text(1.15, 90, 'Flecainide', fontsize=14)
    axs[1].text(1.15, 290, 'Quinine', fontsize=14)

    axs[0].set_ylim(-100, 100)
    axs[1].set_ylim(-100, 300)

    axs[0].set_ylabel(r'$dV/dt_{max}$')

    plt.savefig('./figure-pdfs/f-flec_quin-dvdt-ina.pdf')
    plt.show()


def plot_ap_figure(feature, is_pct_change):
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 4))
    fig.subplots_adjust(.1, .15, .95, .95)

    axs = axs.flatten()

    #is_pct_change = True

    plot_ap_feature(axs[0], 'flecainide', feature, is_pct_change=is_pct_change)
    plot_ap_feature(axs[1], 'quinine', feature, is_pct_change=is_pct_change)

    for vals in [[0, 'A'], [1, 'B']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)
    

    if feature == 'dvdt_max':
        if is_pct_change:
            axs[0].set_ylim(-100, 100)
            axs[1].set_ylim(-100, 300)
            axs[0].text(1.15, 100, 'Flecainide', fontsize=14)
            axs[1].text(1.15, 300, 'Quinine', fontsize=14)
            axs[0].set_ylabel(r'$\Delta$$dV/dt_{max}$ (%)')
        else:
            axs[0].set_ylim(0, 300)
            axs[1].set_ylim(0, 300)
            axs[0].text(1.15, 290, 'Flecainide', fontsize=14)
            axs[1].text(1.15, 290, 'Quinine', fontsize=14)
            axs[0].set_ylabel(r'$\Delta dV/dt_{max}$ (V/s)')


    if feature == 'mdp':
        if is_pct_change:
            axs[0].set_ylim(-15, 30)
            axs[1].set_ylim(-15, 30)
            axs[0].text(1.15, 30, 'Flecainide', fontsize=14)
            axs[1].text(1.15, 30, 'Quinine', fontsize=14)
            axs[0].set_ylabel(r'$\Delta$RMP (%)')

        #axs[0].set_ylim(0, 300)
        #axs[1].set_ylim(0, 300)
        else:
            axs[0].set_ylim(-20, 15)
            axs[1].set_ylim(-20, 15)
            axs[0].text(1.15, 15, 'Flecainide', fontsize=14)
            axs[1].text(1.15, 15, 'Quinine', fontsize=14)
            axs[0].set_ylabel('$\Delta$RMP (mV)')

    if feature == 'apd90':
        if is_pct_change:
            #axs[0].set_ylim(10, 110)
            #axs[1].set_ylim(10, 110)
            axs[0].text(1.15, 100, 'Flecainide', fontsize=14)
            axs[1].text(1.15, 100, 'Quinine', fontsize=14)
            axs[0].set_ylabel(r'$APD_{90}$ (ms)')
        else:
            axs[0].set_ylim(-30, 35)
            axs[1].set_ylim(-30, 35)
            axs[0].text(1.15, 35, 'Flecainide', fontsize=14)
            axs[1].text(1.15, 35, 'Quinine', fontsize=14)
            axs[0].set_ylabel(r'$\Delta APD_{90}$ (ms)')


    plt.savefig(f'./figure-pdfs/fsummary-ap-{feature}.pdf')
    plt.show()


def plot_vc_figure(channel_name, is_change):
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 4))
    fig.subplots_adjust(.1, .15, .95, .95)

    axs = axs.flatten()

    plot_curr(axs[0], 'flecainide', channel_name, is_change=is_change)
    plot_curr(axs[1], 'quinine', channel_name, is_change=is_change)

    for vals in [[0, 'A'], [1, 'B']]:
        axs[vals[0]].set_title(vals[1], y=.94, x=-.15)
        axs[vals[0]].spines['right'].set_visible(False)
        axs[vals[0]].spines['top'].set_visible(False)

    if channel_name == 'Na1_min':
        axs[0].set_ylim(-100, 100)
        axs[1].set_ylim(-100, 300)
        axs[0].text(1.15, 95, 'Flecainide', fontsize=14)
        axs[1].text(1.15, 290, 'Quinine', fontsize=14)
    
    if channel_name == 'Na2_min':
        axs[0].set_ylim(-100, 140)
        axs[1].set_ylim(-100, 140)
        axs[0].text(1.15, 120, 'Flecainide', fontsize=14)
        axs[1].text(1.15, 120, 'Quinine', fontsize=14)

    if is_change:
        axs[0].set_ylabel(f'{channel_name} % Change')
    else:
        axs[0].set_ylabel(f'{channel_name} A/F Change')


    plt.savefig(f'./figure-pdfs/fsummary-vc-{channel_name}.pdf')
    plt.show()


def plot_ap_feature(ax, drug_name, feature_name, is_pct_change=False):
    all_file_names = {}
    all_file_names['dmso'] = get_valid_cells('dmso')
    all_file_names[drug_name] = get_valid_cells(drug_name)

    dmso_dat = []
    drug_dat = []

    for drug, all_f_names in all_file_names.items():
        for f_name in all_f_names:
            ap_summary = pd.read_csv(
                    f'./data/cells/{f_name}/ap_summary_stats.csv')

            # Only select actual APs:
            ap_summary = ap_summary[ap_summary['is_upstroke'] == 1].copy()

            new_list = [f_name, drug, feature_name]

            for comp_name in ['Baseline', '3xEFPC', '10xEFPC',
                                                '20xEFPC', 'Washout']:
                new_list.append(
                        ap_summary[ap_summary['compound'] == comp_name][feature_name].mean())

            if drug == drug_name:
                drug_dat.append(new_list)
            else:
                dmso_dat.append(new_list)

    dmso_df = pd.DataFrame(dmso_dat, columns=['file', 'compound', 'feature', 'baseline', '3x', '10x', '20x', 'wash'])
    drug_df =  pd.DataFrame(drug_dat, columns=['file', 'compound', 'feature', 'baseline', '3x', '10x', '20x', 'wash'])

    if is_pct_change:
        new_drug_cols = [100*(drug_df.iloc[:, 3+i] - drug_df['baseline'])/drug_df['baseline'] for i in range(0, 5)]
        drug_df['baseline'], drug_df['3x'], drug_df['10x'], drug_df['20x'], drug_df['wash'] = new_drug_cols[0], new_drug_cols[1], new_drug_cols[2], new_drug_cols[3], new_drug_cols[4]

        new_dmso_cols = [100*(dmso_df.iloc[:, 3+i] - dmso_df['baseline'])/dmso_df['baseline'] for i in range(0, 5)]
        dmso_df['baseline'], dmso_df['3x'], dmso_df['10x'], dmso_df['20x'], dmso_df['wash'] = new_dmso_cols[0], new_dmso_cols[1], new_dmso_cols[2], new_dmso_cols[3], new_dmso_cols[4]
    else:
        new_drug_cols = [drug_df.iloc[:, 3+i] - drug_df['baseline'] for i in range(0, 5)]
        drug_df['baseline'], drug_df['3x'], drug_df['10x'], drug_df['20x'], drug_df['wash'] = new_drug_cols[0], new_drug_cols[1], new_drug_cols[2], new_drug_cols[3], new_drug_cols[4]

        new_dmso_cols = [dmso_df.iloc[:, 3+i] - dmso_df['baseline'] for i in range(0, 5)]
        dmso_df['baseline'], dmso_df['3x'], dmso_df['10x'], dmso_df['20x'], dmso_df['wash'] = new_dmso_cols[0], new_dmso_cols[1], new_dmso_cols[2], new_dmso_cols[3], new_dmso_cols[4]


    x_arr = np.array([0, 1, 2, 3, 4])
    dmso_means = dmso_df[~dmso_df.isnull().any(axis=1)].iloc[:, 3:].values
    ax.errorbar(x_arr-.05, dmso_means.mean(0), dmso_means.std(0), linestyle='-', c='k', label='DMSO', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08)-.05, y_vals, c='k', alpha=.4) for y_vals in dmso_means]
    drug_means = drug_df[~drug_df.isnull().any(axis=1)].iloc[:, 3:].values
    ax.errorbar(x_arr+0.05, drug_means.mean(0), drug_means.std(0), linestyle='--', c='r', label=drug_name.capitalize(), capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08)+.05, y_vals, c='r', alpha=.4) for y_vals in drug_means]

    #ax.set_ylabel(r'$dV/dt_{max}$ (V/s)')
    #dmso_means = np.array(dmso_means)
    #drug_means = np.array(drug_means)

    max_val = np.concatenate((dmso_means, drug_means)).max()
    min_val = np.concatenate((dmso_means, drug_means)).min()

    y_sig = max_val + np.abs(max_val)*.1
    is_sig = [ttest_ind(dmso_means[:,i], drug_means[:,i]).pvalue < .05 for i in range(0, 4)]
    sig_vals = [ttest_ind(dmso_means[:,i], drug_means[:,i]).pvalue for i in range(0, 4)]

    for i, sig_val in enumerate(sig_vals):
        if sig_val < .001:
            aster = '***'
        elif sig_val < .01:
            aster = '**'
        elif sig_val < .05:
            aster = '*'
        else:
            aster = ''

        ax.text(i-.05, y_sig, aster, fontsize=14)
    
    #[ax.text(i, y_sig, '*', fontsize=14) for i in range(0, 4) if is_sig[i]]


    #ax.legend()
    ax.set_title(drug_name.capitalize())
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(['Baseline', '3xEFPC', '10xEFPC', '20xEFPC', 'Wash'])

    ax.set_ylim(min_val - np.abs(min_val)*.1, max_val + np.abs(max_val*.25))


def plot_curr(ax, drug_name, channel_name, is_change=False):
    all_file_names = {}
    all_file_names['dmso'] = get_valid_cells('dmso')
    all_file_names[drug_name] = get_valid_cells(drug_name)

    dmso_dat = []
    drug_dat = []

    for drug, all_f_names in all_file_names.items():
        for f_name in all_f_names:
            vc_summary = pd.read_csv(
                    f'./data/cells/{f_name}/vc_summary_stats.csv')

            if ('Na1' in channel_name) or ('Na2' in channel_name):
                if (vc_summary[vc_summary['compound'] == 'Baseline']
                                            [channel_name].mean() > -40):
                    continue

            new_list = [f_name, drug, channel_name]

            for comp_name in ['Baseline', '3xEFPC', '10xEFPC',
                                                '20xEFPC', 'Washout']:

                new_list.append(
                        vc_summary[vc_summary['compound'] == comp_name][channel_name].mean())

            if drug == drug_name:
                drug_dat.append(new_list)
            else:
                dmso_dat.append(new_list)

    dmso_dat = pd.DataFrame(dmso_dat, columns=['file', 'compound', 'feature', 'baseline', '3x', '10x', '20x', 'wash'])
    drug_dat =  pd.DataFrame(drug_dat, columns=['file', 'compound', 'feature', 'baseline', '3x', '10x', '20x', 'wash'])

    #import pdb
    #pdb.set_trace()

    if is_change:
        new_drug_cols = [100*(drug_dat.iloc[:, 3+i] - drug_dat['baseline'])/drug_dat['baseline'] for i in range(0, 5)]
        drug_dat['baseline'], drug_dat['3x'], drug_dat['10x'], drug_dat['20x'], drug_dat['wash'] = new_drug_cols[0], new_drug_cols[1], new_drug_cols[2], new_drug_cols[3], new_drug_cols[4]

        new_dmso_cols = [100*(dmso_dat.iloc[:, 3+i] - dmso_dat['baseline'])/dmso_dat['baseline'] for i in range(0, 5)]
        dmso_dat['baseline'], dmso_dat['3x'], dmso_dat['10x'], dmso_dat['20x'], dmso_dat['wash'] = new_dmso_cols[0], new_dmso_cols[1], new_dmso_cols[2], new_dmso_cols[3], new_dmso_cols[4]
    else:
        new_drug_cols = [drug_dat.iloc[:, 3+i] - drug_dat['baseline'] for i in range(0, 5)]
        drug_dat['baseline'], drug_dat['3x'], drug_dat['10x'], drug_dat['20x'], drug_dat['wash'] = new_drug_cols[0], new_drug_cols[1], new_drug_cols[2], new_drug_cols[3], new_drug_cols[4]

        new_dmso_cols = [dmso_dat.iloc[:, 3+i] - dmso_dat['baseline'] for i in range(0, 5)]
        dmso_dat['baseline'], dmso_dat['3x'], dmso_dat['10x'], dmso_dat['20x'], dmso_dat['wash'] = new_dmso_cols[0], new_dmso_cols[1], new_dmso_cols[2], new_dmso_cols[3], new_dmso_cols[4]

    x_arr = np.array([0, 1, 2, 3, 4])



    #FIX THIS DMSO_DIFF ISSUE
    ax.errorbar(x_arr, dmso_dat.iloc[:, 3:].mean(0), dmso_dat.iloc[:, 3:].std(0), linestyle='-', c='k', label='DMSO', capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='k', alpha=.4) for y_vals in dmso_dat.iloc[:, 3:].values]

    x_arr = x_arr + 0.05
    ax.errorbar(x_arr, drug_dat.iloc[:, 3:].mean(0), drug_dat.iloc[:, 3:].std(0), linestyle='--', c='r', label=drug.capitalize(), capsize=3)
    [ax.scatter(x_arr+np.random.uniform(-.08, .08), y_vals, c='r', alpha=.4) for y_vals in drug_dat.iloc[:, 3:].values]

    max_val = np.concatenate((dmso_dat.iloc[:, 3:].values,
                                drug_dat.iloc[:, 3:].values)).max()
    min_val = np.concatenate((dmso_dat.iloc[:, 3:].values,
                                    drug_dat.iloc[:, 3:].values)).min()

    y_sig = max_val + np.abs(max_val)*.1
    is_sig = [ttest_ind(dmso_dat.iloc[:, 3:].values[:, i],
                        drug_dat.iloc[:, 3:].values[:, i]).pvalue < .05
                            for i in range(0, 5)]

    #[ax.text(i, y_sig, '*', fontsize=14) for i in range(0, 5) if is_sig[i]]

    sig_vals = [ttest_ind(dmso_dat.iloc[:, 3:].values[:, i],
                    drug_dat.iloc[:, 3:].values[:, i]).pvalue
                        for i in range(0, 4)]


    for i, sig_val in enumerate(sig_vals):
        if sig_val < .001:
            aster = '***'
        elif sig_val < .01:
            aster = '**'
        elif sig_val < .05:
            aster = '*'
        else:
            aster = ''

        ax.text(i-.05, y_sig, aster, fontsize=14)

    #ax.set_ylabel(r'$I_{Na}$ Change from Baseline (%)')

    #ax.legend()

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(['B', '3x', '10x', '20x', 'W'])
    ax.set_ylim(min_val - np.abs(min_val)*.1, max_val + np.abs(max_val*.25))


# Utility
def moving_average(x, n=10):
    idxs = range(n, len(x), n)
    new_vals = [x[(i-n):i].mean() for i in idxs]
    return np.array(new_vals)


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

    #fig, axs = plt.subplots(3, 1, sharex=True)
    #axs[0].plot(times, stim_region)
    #axs[0].scatter(times_avg, vals)
    #axs[1].plot(times_avg[1:], diff/(times_avg[2]-times_avg[1]))
    #axs[2].plot(times_avg[2:], np.diff(diff))
    #plt.show()

    #import pdb
    #pdb.set_trace()

    return [max_dvdt, is_na]


def get_all_ap_meta(drug_name):
    all_files = get_valid_cells(drug_name)

    meta = {}
    
    for f in all_files:
        #f = 'dmso_220829_001_2' #NO NaV
        ap_dat = pd.read_csv(f'./data/cells/{f}/ap_df.csv')
        ap_meta = pd.read_csv(f'./data/cells/{f}/ap_meta.csv')

        dvdt_is_na = np.array([get_dvdt(k, ap_dat[k]) for k in ap_meta['sweep'].values])
        ap_meta['dvdt_max'] = dvdt_is_na[:, 0]
        ap_meta['is_upstroke'] = dvdt_is_na[:, 1]

        meta[f] = ap_meta

        print(f'\t There are {ap_meta["is_upstroke"].sum()} upstrokes.')
        print(dvdt_is_na[:, 0])

    return meta



def main():
    plot_ap_figure(feature='dvdt_max', is_pct_change=True)
    #plot_ap_figure(feature='mdp',  is_pct_change=False)
    #plot_ap_figure(feature='apd90',  is_pct_change=False)
    #plot_vc_figure('Na1_min', True)#THROW OUT ANY CELLS WITH BASELINE NA >-29 A/F
    #plot_vc_figure('Na2_min', False)#THROW OUT ANY CELLS WITH BASELINE NA >-29 A/F
    #plot_vc_figure('NaL_avg', False)
    #plot_vc_figure('Kr_avg', False)
    #plot_vc_figure('CaL_min', False)
    #plot_vc_figure('to_max', False)
    #plot_vc_figure('K1_avg', False)
    #plot_vc_figure('f1_avg', False)
    #plot_vc_figure('f2_avg', False)
    #plot_vc_figure('Ks1_avg', False)
    #plot_vc_figure('Ks2_avg', False)




if __name__ == '__main__':
    main()

