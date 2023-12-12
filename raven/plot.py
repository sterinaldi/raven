import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import norm
from pathlib import Path
from figaro.plot import plot_median_cr as figaro_plot_median_cr
from figaro import plot_settings

def plot_median_cr(draws, mu, weight, vel_disp, *args, **kwargs):
    single_pdf = lambda x: weight*norm(mu, vel_disp).pdf(x)
    figaro_plot_median_cr(draws, injected = single_pdf, injected_label = '\\sigma = {0:.1f}'.format(vel_disp)+'\ \\mathrm{km/s}', hierarchical = True, *args, **kwargs)
    Path(kwargs.get('out_folder'), 'log_{}.pdf'.format(kwargs.get('name'))).unlink()

def plot_single_fraction(samples, out_folder = '.', name = 'cluster'):
    l_f, m_f, u_f = np.percentile(samples, [16,50,84])
    label = '$w_\\mathrm{single} ='+'{0:.2f}'.format(m_f)+'^{+'+'{0:.2f}'.format(u_f-m_f)+'}_{-'+'{0:.2f}'.format(m_f-l_f)+'}$'
    fig, ax = plt.subplots()
    ax.hist(samples, histtype = 'step', density = True, label = label)
    ax.set_xlabel('$w_\\mathrm{single}$')
    ax.set_ylabel('$p(w_\\mathrm{single})$')
    leg = ax.legend(handlelength=0, handletextpad=0)
    for item in leg.legendHandles:
        item.set_visible(False)
    fig.savefig(Path(out_folder, 'single_fraction_{}.pdf'.format(name)), bbox_inches = 'tight')
    plt.close(fig)

def plot_p_single(probs, stars, threshold, out_folder = '.', name = 'cluster', single = None):
    step = 0.3
    fig, ax = plt.subplots(figsize = (6.4, len(probs)*0.3))
    ax.axvline(threshold, ls = '--', lw = 0.7, c = 'grey', dashes = (5,5), label = '$p_\\mathrm{th} = '+'{}$'.format(threshold))
    if single is None:
        single = (probs > threshold)
    for i, (p, s, single_flag) in enumerate(zip(probs,stars, single)):
        if single_flag:
            marker = '*'
            color  = 'limegreen'
        else:
            marker = 'x'
            color  = 'firebrick'
        if p > 0.5:
            alignment = 'right'
            offset    = -0.02
        else:
            alignment = 'left'
            offset    = 0.02
        ax.scatter([p], [i*step], marker = marker, c = color)
        ax.text(x = p + offset, y = i*step, s = '$\\mathrm{'+'{}'.format(s)+'}$', ha = alignment, va = 'center')
    # Set right limit to 1
    ax.scatter([1], [0], marker = '', c = 'none')
    ax.scatter([0], [0], marker = '', c = 'none')
    ax.set_ylim(-step, (len(probs))*step)
    ax.set_xlabel('$p_\\mathrm{single}$')
    ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    handles, labels = ax.get_legend_handles_labels()
    single_star = Line2D([0],[0], label = '$\\mathrm{Single\ star}$', color = 'limegreen', marker = '*', ls = '')
    binary      = Line2D([0],[0], label = '$\\mathrm{Binary}$', color = 'firebrick', marker = 'x', ls = '')
    handles.extend([single_star, binary])
    ax.legend(loc = 0, handles = handles)
    fig.savefig(Path(out_folder, 'p_single_{}.pdf'.format(name)), bbox_inches = 'tight')
    plt.close(fig)
