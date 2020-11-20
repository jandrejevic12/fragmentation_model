import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
from matplotlib import cm
import seaborn as sns
sns.set_style("darkgrid")
matplotlib.rcParams['grid.linewidth'] = 3

import matplotlib.font_manager as font_manager
font_path = 'Lato-Semibold.ttf'
prop = font_manager.FontProperties(fname=font_path)

cmap = matplotlib.cm.get_cmap("plasma")
cnorm = matplotlib.colors.Normalize(vmin=0, vmax=1)

def get_cmap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def fig_with_cbar(figsize, ncols):
    fig = plt.figure(figsize=figsize)
    nrows = 1
    gs = GridSpec(nrows,ncols+1,width_ratios=[25]*ncols+[1])
    axes = [fig.add_subplot(gs[i]) for i in range(nrows*ncols)]
    cbax = fig.add_subplot(gs[0,-1])
    return fig, axes, cbax

def set_axis_labels(ax, xlabel, ylabel, size, xticks=False, yticks=False):
    ax.tick_params(axis='both', which='major', pad=12)
    #if xticks and xlabel:
    if xticks:
        ax.set_xticks(xticks)
    #if yticks and ylabel:
    if yticks:    
        ax.set_yticks(yticks)
    if xlabel:
        prop.set_size(size)
    else:
        prop.set_size(0)
    for lab in ax.get_xticklabels():
        lab.set_fontproperties(prop)
    ax.set_xlabel(xlabel, fontproperties=prop)
    if ylabel:
        prop.set_size(size)
    else:
        prop.set_size(0)
    for lab in ax.get_yticklabels():
        lab.set_fontproperties(prop)
    ax.set_ylabel(ylabel, fontproperties=prop)

def set_legend(ax, loc, size, use_prop=True):
    prop.set_size(size)
    if use_prop:
        ax.legend(loc=loc, prop=prop, frameon=False)
    else:
        ax.legend(loc=loc, fontsize=size, frameon=False)

def set_text(ax, x, y, text, size):
    prop.set_size(size)
    ax.text(x, y, text, fontproperties=prop)

def set_cbar(fig, cbax, cmap, cnorm, label, ticks, size):
    prop.set_size(size)
    sm = cm.ScalarMappable(cmap=cmap, norm=cnorm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbax)
    cbar.set_label(label, fontproperties=prop)
    cbar.ax.get_yaxis().set_ticks(ticks)
    cbar.ax.set_yticklabels([str(t) for t in ticks])
    for lab in cbar.ax.get_yticklabels():
        lab.set_fontproperties(prop)

