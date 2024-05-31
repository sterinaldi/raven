import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rcParams
from scipy.stats import norm
from pathlib import Path
from figaro.plot import plot_median_cr as figaro_plot_median_cr
from figaro import plot_settings
rcParams["axes.grid"] = False

def plot_median_cr(draws, mu, weight, vel_disp, *args, **kwargs):
    single_pdf = lambda x: weight*norm(mu, vel_disp).pdf(x)
    figaro_plot_median_cr(draws, injected = single_pdf, injected_label = '\\sigma_\\mathrm{V} = '+'{0:.1f}'.format(vel_disp)+'\ \\mathrm{km/s}', hierarchical = True, *args, **kwargs)
    # Removing unwanted FIGARO plot
    Path(kwargs.get('out_folder'), 'log_observed_{}.pdf'.format(kwargs.get('name'))).unlink()
    # Removing observed from file names
    Path(kwargs.get('out_folder'), 'observed_{}.pdf'.format(kwargs.get('name'))).rename(Path(kwargs.get('out_folder'), '{}.pdf'.format(kwargs.get('name'))))
    Path(kwargs.get('out_folder'), 'prob_observed_{}.txt'.format(kwargs.get('name'))).rename(Path(kwargs.get('out_folder'), 'prob_{}.txt'.format(kwargs.get('name'))))

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

def plot_p_single(probs, stars, out_folder = '.', name = 'cluster', single = None):
    step = 0.4
    fig, ax = plt.subplots(figsize = (6.4, len(probs)*step))
    # Set right limit to 1
    ax.scatter([1], [-step], marker = '', c = 'none')
    ax.scatter([0], [-step], marker = '', c = 'none')
    xmin, xmax = ax.get_xlim()
    # Delimiters
    ax.axvline(0.5, ls = '--', lw = 0.7, c = 'grey', dashes = (5,5), alpha = 0.7)
    ax.axvline(0.1, ls = '--', lw = 0.7, c = 'grey', dashes = (5,5), alpha = 0.7)
    ax.axvline(0.9, ls = '--', lw = 0.7, c = 'grey', dashes = (5,5), alpha = 0.7)
    ax.text(x = xmin +(0.1-xmin)/2., y = -0.3, s = '$\\mathrm{Binary}$', ha = 'center', va = 'center', color = 'grey')
    ax.text(x = 0.3, y = -0.3, s = '$\\mathrm{Potential}\\ \mathrm{binary}$', ha = 'center', va = 'center', color = 'grey')
    ax.text(x = 0.7, y = -0.3, s = '$\\mathrm{Potential}\\ \mathrm{single}$', ha = 'center', va = 'center', color = 'grey')
    ax.text(x = xmax - (xmax-0.9)/2., y = -0.3, s = '$\\mathrm{Single}$', ha = 'center', va = 'center', color = 'grey')
    if single is None:
        single_flag = {s: None for s in stars}
    else:
        single_flag = single
    NA_legend     = False
    binary_legend = False
    single_legend = False
    sana_flag = False
    for i, (p, s) in enumerate(zip(probs, stars)):
        ll, l, m, u, uu = np.percentile(p, [5,16,50,84,95])
        if single is not None:
            sana_flag = True
            if single_flag[s] is not None:
                if single_flag[s]:
                    single_legend = True
                    marker    = '*'
                    color     = 'limegreen'
                    facecolor = color
                    edgecolor = color
                else:
                    binary_legend = True
                    marker    = 'x'
                    color     = 'firebrick'
                    facecolor = color
                    edgecolor = None
            else:
                NA_legend = True
                marker    = 'o'
                color     = 'steelblue'
                facecolor = 'w'
                edgecolor = color
        else:
            marker    = 'o'
            color     = 'steelblue'
            facecolor = 'w'
            edgecolor = color
        if m < 0.1:
            textpos = 0
            ha = 'left'
        elif m > 0.9:
            textpos = 1
            ha = 'right'
        else:
            textpos = m
            ha = 'center'
        ax.errorbar([m], [i*step], xerr = [[m-ll],[uu-m]], marker = '', color = color, alpha = 0.5, lw = 1, capsize = 2)
        ax.errorbar([m], [i*step], xerr = [[m-l],[u-m]], marker = '', color = color, alpha = 0.75, lw = 1, capsize = 2)
        ax.scatter([m], [i*step], marker = marker, edgecolors = edgecolor, facecolors = facecolor, zorder = 3*len(stars)+10, lw = 1)
        ax.text(x = textpos, y = i*step+0.175, s = '$\\mathrm{'+'{}'.format(s)+'}$', ha = ha, va = 'center')
    # Set limits
    ax.set_ylim(-0.6, (len(probs))*step)
    ax.set_xlabel('$p_\\mathrm{single}$')
    ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    if sana_flag:
        handles, labels = ax.get_legend_handles_labels()
        if single_legend:
            single_star = Line2D([0],[0], label = '$\\mathrm{Single\ star}$', color = 'limegreen', marker = '*', ls = '')
            handles.extend([single_star])
        if binary_legend:
            binary      = Line2D([0],[0], label = '$\\mathrm{Binary}$', color = 'firebrick', marker = 'x', ls = '')
            handles.extend([binary])
        if NA_legend:
            NA_handle = Line2D([0],[0], label = '$\\mathrm{N/A\ (single\ epoch)}$', color = 'steelblue', marker = 'o', ls = '', markerfacecolor = 'white')
            handles.extend([NA_handle])
        if (single_legend or binary_legend):
            ax.legend(loc = 0, handles = handles, fontsize = 10, title = '$\\mathrm{Method\ from\ Sana\ et\ al.\ (2013):}$')
    fig.savefig(Path(out_folder, 'p_single_{}.pdf'.format(name)), bbox_inches = 'tight')
    plt.close(fig)
