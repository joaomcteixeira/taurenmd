"""Plot a single parameter."""
import itertools
import numpy as np

from matplotlib import pyplot as plt

from bioplottemplates import log
from bioplottemplates.libs import libmsg, libutil
from bioplottemplates.logger import S, T

def plot(
        x_labels,
        y_data,
        title=None,
        xlabel=None,
        ylabel=None,
        series_labels=None,
        legend=True,
        legend_fs=6,
        legend_loc=4,
        numeric_x_labels=False,
        colors=('b', 'g', 'r', 'c', 'm', 'y', 'k'),
        alpha=0.7,
        grid=True,
        grid_color="lightgrey",
        grid_ls="-",
        grid_lw=1,
        grid_alpha=0.5,
        figsize=(10, 6),
        filename='plot_param.pdf',
        **kwargs,
        ):
    """
    """
    log.info(T('Plotting Labeled Dots'))
    
    # prepares data
    if isinstance(y_data, (list, np.ndarray)):
        if not isinstance(y_data[0], (list, np.ndarray)):
            y_data = [y_data]
    else:
        raise ValueError('y_data must be list or np.ndarray')
    
    if not isinstance(series_labels, list):
        series_labels = [series_labels]
    elif series_labels is None:
        series_labels = [None] * len(y_data)

    plot_colors = itertools.cycle(libutil.make_list(colors))
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
    if xlabel:
        plt.tight_layout(rect=[0.05, 0.15, 0.995, 0.985])
    else:
        plt.tight_layout(rect=[0.05, 0.10, 0.995, 0.985])

    for dataset in y_data:
        fig.suptitle(
            title,
            x=0.5,
            y=0.990,
            va="top",
            ha="center",
            )
    
    #ax.margins(x=1)
    for i, yy in enumerate(y_data):

        ax.scatter(
            range(int(numeric_x_labels), int(numeric_x_labels) + len(yy)),
            yy,
            label=series_labels[i],
            color=next(plot_colors),
            alpha=alpha,
            zorder=10)
    
    
    ax.set_xlabel(xlabel, weight='bold')
    ax.set_ylabel(ylabel, weight='bold')
    
    #ax.set_xlim(x_data[0], x_data[-1])
    ax.set_ylim(0)
    
    if numeric_x_labels:
        ax.set_xlim(numeric_x_labels, len(x_labels))
    else:
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=90)

    if grid:
        ax.grid(
            which='major',
            color=grid_color,
            linestyle=grid_ls,
            linewidth=grid_lw,
            alpha=grid_alpha,
            zorder=1,
            )
    
    if legend:
        ax.legend(
            fontsize=legend_fs,
            loc=legend_loc,
            )

    fig.savefig(filename)
    
