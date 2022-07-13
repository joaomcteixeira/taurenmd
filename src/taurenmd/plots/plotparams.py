"""Plot a single parameter."""
import itertools

import numpy as np
from matplotlib import pyplot as plt

from taurenmd.libs import libutil


def plot(
        x_data,
        y_data,
        *,
        labels="No label provided",
        title=None,
        xlabel=None,
        ylabel=None,
        xlabel_fs=8,
        ylabel_fs=8,
        xticks=None,
        yticks=None,
        xticks_labels=None,
        yticks_labels=None,
        colors=('b', 'g', 'r', 'c', 'm', 'y', 'k'),
        alpha=0.7,
        xmax=None,
        xmin=None,
        ymax=None,
        ymin=None,
        hline=None,
        hline_color="black",
        hline_lw=1,
        xticks_fs=10,
        yticks_fs=10,
        grid=True,
        grid_color="lightgrey",
        grid_ls="-",
        grid_lw=1,
        grid_alpha=0.5,
        legend=True,
        legend_fs=8,
        legend_loc=0,
        vert_lines=None,
        figsize=(8, 5),
        filename='plot_param.pdf',
        dpi=150,
        ):
    """
    Plot a parameter.

    Bellow parameters concern data representation and are considered
    of highest importance because their incorrect use can mislead
    data analysis and consequent conclusions.

    Plot style parameters concernning only plot style, i.e., colors,
    shapes, fonts, etc... and which do not distort the actual data,
    are not listed in the paremeter list bellow. We hope these
    parameter names are self-explanatory and are listed in the function
    definition.

    Parameters
    ----------
    x_data : interable of numbers
        Container of the X axis data. Should be accepted
        by matplotlib.

    y_data : np.ndarray, shape=(M, len(x_data))
        Container of the Y axis data.
        Where M is the number of series.

    labels : str, optional
        The label to represent in plot legend.
        If a list of series is provided, a list of labels
        can be provided as well.
        Defauts to: "no labels provided".

    filename : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to rmsd_individual_chains_one_subplot.pdf.
        You can change the file type by specifying its extention in
        the file name.

    fig_size : tuple of float or int
        The size ratio of the subplot in the figure.
    """
    y_data = np.array(y_data)

    if y_data.ndim == 1:
        y_data = y_data[np.newaxis, :]

    plot_labels = libutil.make_list(labels)
    plot_colors = libutil.make_list(colors)
    plot_colors = itertools.cycle(plot_colors)

    assert len(x_data) == y_data.shape[1], \
        '{} vs {}'.format(len(x_data), y_data.shape[1])
    assert y_data.shape[0] == len(plot_labels), (
        '{} vs {}'.format(y_data.shape[0], len(plot_labels)))

    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=figsize,
        constrained_layout=False,
        )
    plt.tight_layout(rect=[0.05, 0.05, 0.995, 1])

    ax.set_title(
        title,
        x=0.5,
        y=1,
        va="bottom",
        ha="center",
        fontweight='bold',
        )

    for yy, ll in zip(y_data, plot_labels):
        ax.plot(
            x_data,
            yy,
            label=ll,
            color=next(plot_colors),
            alpha=alpha,
            zorder=2,
            )

    ax.set_xlabel(xlabel, weight='bold', fontsize=xlabel_fs)
    ax.set_ylabel(ylabel, weight='bold', fontsize=ylabel_fs)

    # setting axes scales

    xmin = xmin if xmin is not None else x_data[0]
    xmax = xmax if xmax is not None else x_data[-1]
    ymin = ymin if ymin is not None else np.array(y_data).min()
    ymax = ymax if ymax is not None else np.array(y_data).max()

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    if xticks:
        ax.set_xticks(ticks=xticks)

    if xticks_labels:
        ax.set_xticklabels(xticks_labels)

    if yticks:
        ax.set_yticks(ticks=yticks)

    if yticks_labels:
        ax.set_yticklabels(yticks_labels)

    ax.tick_params(axis="x", labelsize=xticks_fs)
    ax.tick_params(axis="y", labelsize=yticks_fs)

    if grid:
        ax.grid(
            color=grid_color,
            linestyle=grid_ls,
            linewidth=grid_lw,
            alpha=grid_alpha,
            zorder=1,
            )

    if hline is not None:
        ax.axhline(hline, zorder=1, color=hline_color, lw=hline_lw)

    if isinstance(vert_lines, (list, tuple)):
        for line in vert_lines:
            ax.axvline(x=float(line), color='k')

    if legend:
        ax.legend(
            fontsize=legend_fs,
            loc=legend_loc,
            )
    plt.subplots_adjust(top=.9)
    fig.savefig(filename, dpi=dpi)

    plt.close("all")

    return
