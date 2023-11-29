import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from utility import get_valid_cells
from scipy.stats import ttest_ind

plt.rcParams['lines.linewidth'] = .9
plt.rcParams['lines.markersize'] = 4
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rc('legend', fontsize = 8)
matplotlib.rcParams['pdf.fonttype'] = 42

def plot_figure():
    # quinine_221206_003_1
    # quinine_221206_009_1
    # quinine_221206_009_2
    # quinine_230120_007_3
    # quinine_230120_010_2
    # quinine_230403_011_1
    # quinine_230403_011_3
    fig = plt.figure(figsize=(6.5, 6))
    fig.subplots_adjust(.1, .1, .95, .93)
    grid = fig.add_gridspec(4, 4, hspace=.6, wspace=.5)

    sub = grid[0, :]
    sub_1 = sub.subgridspec(1, 4, hspace=.3)
    axs_1 = [fig.add_subplot(sub_1[i]) for i in range(0, 4)]

    sub = grid[1, :]
    sub_2 = sub.subgridspec(1, 4, hspace=.3)
    axs_2 = [fig.add_subplot(sub_2[i]) for i in range(0, 4)]

    sub = grid[2:, 0:2]
    sub_flec = sub.subgridspec(1, 1)
    ax_flec = fig.add_subplot(sub_flec[0])

    sub = grid[2:, 2:]
    sub_quin = sub.subgridspec(1, 1)
    ax_quin = fig.add_subplot(sub_quin[0])

    plot_ap_dat(axs_1, 'flecainide_221117_001_3', False)
    plot_ap_dat(axs_2, 'quinine_230403_011_3')

    #plot_upstroke_change(ax_flec, 'flecainide')
    #plot_upstroke_change(ax_quin, 'quinine')
    plot_ap_feature(ax_flec, 'flecainide', 'dvdt_max', is_pct_change=True)
    plot_ap_feature(ax_quin, 'quinine', 'dvdt_max', is_pct_change=True)

    letters = ['C', 'D']
    for i, ax in enumerate([ax_flec, ax_quin]):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_title(letters[i], y=.94, x=-.1)
    
    ax_flec.set_xlabel('Flecainide Conc.')
    ax_quin.set_xlabel('Quinine Conc.')
    
    axs_1[0].text(-7, 50, 'A')
    axs_1[0].text(-2, 30, 'Flec')

    axs_2[0].text(-7, 50, 'B')
    axs_2[0].text(-2, 30, 'Quin')
    
    ax_flec.set_ylabel(r'Change in $dV/dt_{max}$ (%)')


    plt.savefig('./figure-pdfs/f-ap-upstroke-drug.pdf')
    plt.show()


def plot_ap_dat(axs, cell_folder, with_xlabel=True):
    ap_dat = pd.read_csv(f'data/cells/{cell_folder}/ap_df.csv')
    ap_meta = pd.read_csv(f'data/cells/{cell_folder}/ap_meta.csv')

    trial_nums = np.sort(np.unique([int(trial.split('_')[1]) for trial in ap_dat.keys()]))


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


    for i, ax in enumerate(axs):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        if with_xlabel:
            ax.set_xlabel('Time (ms)')

        #ax.arrow(2.5, 20, 5, -30, head_width=5, head_length=10, fc='k')
        if i == 3:
            ax.annotate("", xy=(1, 20), xytext=(7, -20), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='grey'))
        else:
            ax.annotate("", xy=(7, -20), xytext=(1, 20), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='red'))

    if not with_xlabel:
        axs[0].set_title('3xEFPC')
        axs[1].set_title('10xEFPC')
        axs[2].set_title('20xEFPC')
        axs[3].set_title('DMSO wash')
    axs[0].set_ylabel('mV')
    #axs[2].set_ylabel('mV')

    for ax in axs:
        ax.set_xlim(-3, 12)


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


def main():
    plot_figure()

if __name__ == '__main__':
    main()