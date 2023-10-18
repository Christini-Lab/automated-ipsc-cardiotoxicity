import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
    # flecainide_220901_009_3
    # flecainide_221117_001_2
    # 
    plot_ap_dat('flecainide_221117_001_3')


def plot_ap_dat(cell_folder):
    ap_dat = pd.read_csv(f'data/cells/{cell_folder}/ap_df.csv')
    ap_meta = pd.read_csv(f'data/cells/{cell_folder}/ap_meta.csv')

    trial_nums = np.sort(np.unique([int(trial.split('_')[1]) for trial in ap_dat.keys()]))

    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6), sharey=True, sharex=True)
    fig.subplots_adjust(wspace=.1, hspace=.2)
    axs = axs.flatten()

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

        ax.set_xlabel('Time (ms)')
        #ax.arrow(2.5, 20, 5, -30, head_width=5, head_length=10, fc='k')
        if i == 3:
            ax.annotate("", xy=(3, 20), xytext=(7, -20), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='grey'))
        else:
            ax.annotate("", xy=(7, -20), xytext=(3, 20), arrowprops=dict(headwidth=7, headlength=7, width=0.1, color='red'))

    axs[0].set_title('3xEFPC')
    axs[1].set_title('10xEFPC')
    axs[2].set_title('20xEFPC', y=.9)
    axs[3].set_title('DMSO wash', y=.9)
    axs[0].set_ylabel('Voltage (mV)')
    axs[2].set_ylabel('Voltage (mV)')

    axs[0].set_xlim(-3, 12)

    plt.savefig('./figure-pdfs/f-ap-upstroke-ex.pdf')
    plt.show()



def create_figure_with_four_axes():
    # Create a figure and a 2x2 grid of axes
    fig, axs = plt.subplots(2, 2, figsize=(8, 6))

    # Remove the border on the top and right axes
    axs[0, 0].spines['top'].set_visible(False)
    axs[0, 0].spines['right'].set_visible(False)

    axs[0, 1].spines['top'].set_visible(False)
    axs[0, 1].spines['right'].set_visible(False)

    axs[1, 0].spines['top'].set_visible(False)
    axs[1, 0].spines['right'].set_visible(False)

    axs[1, 1].spines['top'].set_visible(False)
    axs[1, 1].spines['right'].set_visible(False)

    return fig, axs


def main():
    plot_figure()

if __name__ == '__main__':
    main()