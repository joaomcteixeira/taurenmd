"""
Tauren-MD Plotting Functions.
"""
# Copyright © 2018-2019 Tauren-MD Project
#
# Tauren-MD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Tauren-MD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Tauren-MD. If not, see <http://www.gnu.org/licenses/>.
#
# Contributors to this file:
# - João M.C. Teixeira (https://github.com/joaomcteixeira)
import itertools as it
import numpy as np
from functools import wraps

from matplotlib import pyplot as plt
from matplotlib import colors as mcolors

from tauren import logger

log = logger.get_log(__name__)

_msg_fig_saved = "* Saved figure: {}"


def _check_data(func):
    
    @wraps(func)
    def wrapper(*args, **kwargs):
        
        if not(isinstance(args[1], np.ndarray)) \
                or not(isinstance(args[0], np.ndarray)):
            log.info(
                "* Data Type Error * "
                f"{func.__name__} receives np.ndarrays as"
                " second data parameter."
                f" {type(args[1])} given."
                )
            return
        
        if args[0].size == 0 or args[1].size == 0:
            log.info(
                "* Data Error *"
                f" No data given in {func.__name__} function."
                )
            return
        
        result = func(*args, **kwargs)
        
        return result
    
    return wrapper
        

def _calc_fig_size(
        nsubplots,
        *,
        ncols=1,
        irow=4.8,
        icol=6.4,
        ):
    """
    Calculates the figure dimensions (size in inches).
    
    Figure dimensions are calculated based on the following parameters:
    
    Parameters
    ----------
    nsubplots : int
        The total number of subplots
    
    ncols : int, optional, defaults 1
        The desired number of columns
    
    irow : float, optional, defaults 4.8
        Number os inches per row.
    
    icol : float, optional, defaults 6.4
        Number of inches per column.
    
    
    Returns
    -------
    A tuple (width, height) in inches.
    """
    
    # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.figure.html
    
    width = ncols * icol
    height = (nsubplots // ncols) * irow
    
    fs = (width, height)
    
    log.debug(f"returning: {fs}")
    
    return fs


@_check_data
def rmsd_chain_per_subplot(
        x_data,
        y_data,
        *,
        labels="No labels provided",
        suptitle="RMSDs per chain",
        x_label="Frame Number",
        y_label="RMSDs",
        colors=list(mcolors.BASE_COLORS.keys()),
        alpha=0.7,
        grid=True,
        grid_color="lightgrey",
        grid_ls="-",
        grid_lw=1,
        grid_alpha=0.5,
        legend=True,
        legend_fs=6,
        legend_loc=4,
        fig_name="rmsd_chain_per_subplot.pdf",
        **kwargs
        ):
    """
    Plots a single plot with the combined RMSD.
    
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
    x_data : interable
        Container of the X axis data. Should be accepted
        by matplotlib.
    
    y_data : np.ndarray, shape=(N, M)
        Container of the Y axis data.
        Where N is the number of chains (data series), and
        M the data for each series (RMSDs of a given chain)
    
    labels : str, optional
        The chain labels to represent in plot legend.
        Defauts to: "no labels provided".
    
    fig_name : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to rmsd_individual_chains_one_subplot.pdf.
        You can change the file type by specifying its extention in
        the file name.
    """
    log.info("* Plotting RMSDs per chain...")
    
    figsize = _calc_fig_size(
        y_data.shape[0],
        ncols=1,
        irow=3,
        )
    
    fig, axs = plt.subplots(
        nrows=y_data.shape[0],
        ncols=1,
        figsize=figsize,
        sharex=True,
        )
    
    fig.suptitle(
        suptitle,
        x=0.5,
        y=0.990,
        va="top",
        ha="center",
        )
    
    try:
        axs = axs.ravel()
    
    except AttributeError:
        axs = np.array([axs])
        # pass
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    fig.subplots_adjust(hspace=0)
    
    colors_it = it.cycle(colors)
    
    max_rmsd = 0
    for ii, chain_rmsds in enumerate(y_data):
        
        axs[ii].plot(
            x_data,
            chain_rmsds,
            label=labels[ii],
            color=next(colors_it),
            alpha=alpha,
            )
        
        axs[ii].set_xlim(x_data[0], x_data[-1])
        
        axs[ii].set_ylabel(y_label, weight='bold')
        
        if grid:
            axs[ii].grid(
                color=grid_color,
                linestyle=grid_ls,
                linewidth=grid_lw,
                alpha=grid_alpha,
                )
    
        if legend:
            axs[ii].legend(
                fontsize=legend_fs,
                loc=legend_loc,
                )
        
        if chain_rmsds.max() > max_rmsd:
            max_rmsd = chain_rmsds.max()
        
    else:
        axs[ii].set_xlabel(x_label, weight='bold')
    
    log.debug(f"<max_rmsds>: {max_rmsd}")
    
    for axis in axs:
        
        axis.set_ylim(0, max_rmsd * 1.1)
        all_ticks = axis.get_yticks()
        axis.set_yticks(all_ticks[1:-1])
    
    else:
        axis.set_yticks(all_ticks[:-1])
    
    fig.savefig(fig_name)
    log.info(_msg_fig_saved.format(fig_name))
    
    plt.close("all")
    
    return


@_check_data
def rmsd_combined_chains(
        x_data,
        y_data,
        *,
        label="No labels provided",
        suptitle="Combined Chain RMSDs",
        x_label="Frame Number",
        y_label="RMSDs",
        color='blue',
        alpha=0.7,
        grid=True,
        grid_color="lightgrey",
        grid_ls="-",
        grid_lw=1,
        grid_alpha=0.5,
        legend=True,
        legend_fs=6,
        legend_loc=4,
        fig_name='plot_rmsd_combined.pdf',
        **kwargs
        ):
    """
    Plots a single plot with the combined RMSD.
    
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
    
    y_data : np.ndarray, shape=(M,)
        Container of the Y axis data.
        Where M is the RMSDs data for the combined chains.
    
    label : str, optional
        The label to represent in plot legend.
        Defauts to: "no labels provided".
    
    fig_name : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to rmsd_individual_chains_one_subplot.pdf.
        You can change the file type by specifying its extention in
        the file name.
    """
    log.info("* Plotting combined Chain RMSDs...")
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    
    fig.suptitle(
        suptitle,
        x=0.5,
        y=0.990,
        va="top",
        ha="center",
        )
    
    ax.plot(
        x_data,
        y_data,
        label=label,
        color=color,
        alpha=alpha,
        )
    
    ax.set_xlabel(x_label, weight='bold')
    ax.set_ylabel(y_label, weight='bold')
    
    ax.set_xlim(x_data[0], x_data[-1])
    ax.set_ylim(0)
    
    if grid:
        ax.grid(
            color=grid_color,
            linestyle=grid_ls,
            linewidth=grid_lw,
            alpha=grid_alpha,
            )
    
    if legend:
        ax.legend(
            fontsize=legend_fs,
            loc=legend_loc,
            )
    
    fig.savefig(fig_name)
    log.info(_msg_fig_saved.format(fig_name))
        
    plt.close("all")
    
    return


@_check_data
def rmsd_individual_chains_one_subplot(
        x_data,
        y_data,
        *,
        labels="No labels provided",
        suptitle="Chains' RMSDs",
        x_label="Frame Number",
        y_label="RMSDs",
        colors=list(mcolors.BASE_COLORS.keys()),
        alpha=0.7,
        grid=True,
        grid_color="lightgrey",
        grid_ls="-",
        grid_lw=1,
        grid_alpha=0.5,
        legend=True,
        legend_fs=6,
        legend_loc=4,
        fig_name="rmsd_individual_chains_one_subplot.pdf",
        **kwargs
        ):
    """
    Plots a single plot with the combined RMSD.
    
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
    x_data : interable
        Contains X axis data.
        Should be accepted by matplotlib.
    
    y_data : np.ndarray, shape=(N, M)
        Container of the Y axis data.
        Where N is the number of chains (data series), and
        M the data for each series (RMSDs of a given chain)
    
    labels : str, optional
        The chain labels to represent in plot legend.
        Defauts to: "no labels provided".
    
    fig_name : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to rmsd_individual_chains_one_subplot.pdf.
        You can change the file type by specifying its extention in
        the file name.
    """
    
    log.info("* Plotting chains RMSDs single subplot...")
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    fig.suptitle(
        suptitle,
        x=0.5,
        y=0.990,
        va="top",
        ha="center",
        )
    
    colors_it = it.cycle(colors)
    max_rmsds = 0
    
    for ii, chain_rmsd_data in enumerate(y_data):
        
        ax.plot(
            x_data,
            chain_rmsd_data,
            label=labels[ii],
            color=next(colors_it),
            alpha=alpha,
            )
        
        if chain_rmsd_data.max() > max_rmsds:
            max_rmsds = chain_rmsd_data.max()
    
    ax.set_xlim(x_data[0], x_data[-1])
    ax.set_ylim(0, max_rmsds * 1.1)
    
    ax.set_xlabel(x_label, weight='bold')
    ax.set_ylabel(y_label, weight='bold')
    
    if grid:
        ax.grid(
            color=grid_color,
            linestyle=grid_ls,
            linewidth=grid_lw,
            alpha=grid_alpha,
            )
    
    if legend:
        ax.legend(
            fontsize=legend_fs,
            loc=legend_loc,
            )
    
    fig.savefig(fig_name)
    log.info(_msg_fig_saved.format(fig_name))
    
    plt.close("all")
    
    return
