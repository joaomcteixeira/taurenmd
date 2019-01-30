"""
Fuctions to plot RMSDs.

Copyright © 2018-2019 Tauren-MD Project

Contributors to this file:
- João M.C. Teixeira (https://github.com/joaomcteixeira)

Tauren-MD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Tauren-MD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tauren-MD. If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import numpy as np
from matplotlib import pyplot as plt

from tauren import tlog
from tauren._core import validators
from tauren.plot import _commons as pltcommons
from tauren.calculate import fromtraj

log = tlog.get_log(__name__)


@validators.validate_trajectory
def plot_rmsd_combined(
        traj,
        *args,
        color='blue',
        chains="all",
        fig_name='plot_rmsd_combined.pdf',
        **kwargs,
        ):
    """
    Plots a single plot with the combined RMSD for all chains.
    
    Parameters:
    
        - color (str, opt): plot color
        
        - chains (str, opt, def: "all"): the chains to consider for
            plotting. Accepts integer slicer: 1:, :50, 2:30:2, etc...
        
        - fig_name (str, opt, def: "RMSDS_plot_combined.pdf"): the
            name of the output figure.
    """
    
    log.info("* Plotting combined RMSDs...")
    
    validators.validate_args(
        (
            (color, str),
            (chains, str),
            (fig_name, str),
            ),
        func_name=plot_combined.__name__,
        )
    
    if chains != "all":
        chain_list = pltcommons._get_chain_list(traj, chains)
        selector_traj_str = [f"chainid {chain}" for chain in chain_list]
        selector = " or ".join(selector_traj_str)
        slicer = traj.topology.select(selector)
        try:
            sliced_traj = traj.atom_slice(slicer, inplace=False)
        except IndexError:
            log.expection("Could not slice traj")
            sys.exit(1)
        
        _, rmsds = fromtraj.calc_rmsds(sliced_traj)
    
    else:
        _, rmsds = fromtraj.calc_rmsds(traj)
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    fig.suptitle("Combined RMSDs", x=0.5, y=0.990, va="top", ha="center")
    
    x_range = np.arange(0, traj.n_frames)
    
    ax.plot(
        x_range,
        rmsds,
        label=chains,
        color=color,
        alpha=0.7
        )
    
    ax.set_xlabel("Frame Number", weight='bold')
    ax.set_ylabel("RMSDs", weight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=1, alpha=0.7)
    
    ax.set_xlim(x_range[0], x_range[-1])
    ax.set_ylim(0)
    
    ax.legend(fontsize=6, loc=4)
    
    fig.savefig(fig_name)
    log.info(f"    saved figure '{fig_name}'")
    
    return (traj, )


@validators.validate_trajectory
def plot_rmsd_chain_per_subplot(
        traj,
        *args,
        chains="all",
        colors="None",
        fig_name="plot_rmsd_chain_per_subplot.pdf",
        **kwargs,
        ):
    """
    Plots RMSD for individual chains.
    
    Draws a figure with one subplot for each chain.
    
    Parameters:
    
        - chains (str, opt, def: "all"): the chains to consider for
            plotting. Accepts integer slicer: 1:, :50, 2:30:2, etc...
        
        - colors (str of comma separated values): "blue,red,green".
        
        - fig_name (str, opt, def: "RMSDS_plot_combined.pdf"): the
            name of the output figure.
    """
    
    log.info("* Plotting RMSDs per chain...")
    
    validators.validate_args(
        (
            (chains, str),
            (colors, str),
            (fig_name, str),
            ),
        func_name=plot_chain_per_subplot.__name__,
        )
    
    chain_list = pltcommons._get_chain_list(traj, chains)
    colors = pltcommons._get_colors(colors)
    figsize = pltcommons._fig_size(len(chain_list), ncols=1, irow=3)
    
    fig, ax = plt.subplots(
        nrows=len(chain_list),
        ncols=1,
        figsize=figsize,
        sharex=True,
        )
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    fig.suptitle("RMSDs per chain", x=0.5, y=0.990, va="top", ha="center")
    
    ax = ax.ravel()
    
    # https://matplotlib.org/gallery/subplots_axes_and_figures/ganged_plots.html#sphx-glr-gallery-subplots-axes-and-figures-ganged-plots-py
    fig.subplots_adjust(hspace=0)
    
    x_range = np.arange(0, traj.n_frames)
    
    max_rmsd = 0
    for ii, chain in enumerate(chain_list):
        
        selector = f"chainid {chain}"
        slicer = traj.topology.select(selector)
        log.debug(
            f"performing on: selector {selector} "
            f"with slicer len {len(slicer)}"
            )
        
        if slicer.size == 0:
            errmsg = f"chainid '{chain}' does not exists. Ignoring..."
            log.info(errmsg)
            ax[ii].text(
                0.5,
                0.5,
                errmsg,
                fontsize=10,
                va='center',
                ha='center',
                transform=ax[ii].transAxes
                )
            rmsds = np.zeros(traj.n_frames)
            continue
        
        else:
            chain_traj = traj.atom_slice(slicer, inplace=False)
            _, rmsds = fromtraj.calc_rmsds(chain_traj)
        
        ax[ii].plot(
            x_range,
            rmsds,
            label=chain,
            color=next(colors),
            alpha=0.7,
            )
        
        ax[ii].set_xlim(x_range[0], x_range[-1])
        
        ax[ii].set_ylabel("RMSDs", weight='bold')
        
        ax[ii].grid(color='lightgrey', linestyle='-', linewidth=1, alpha=0.7)
        ax[ii].legend(fontsize=6, loc=4)
        
        if rmsds.max() > max_rmsd:
            max_rmsd = rmsds.max()
        
    else:
        ax[ii].set_xlabel("Frame Number", weight='bold')
    
    log.debug(f"<max_rmsds>: {max_rmsd}")
    
    for i, z in enumerate(ax):
        
        ax[i].set_ylim(0, max_rmsd + max_rmsd * 0.1)
        all_ticks = ax[i].get_yticks()
        ax[i].set_yticks(all_ticks[1:-1])
    
    else:
        ax[i].set_yticks(all_ticks[:-1])
    
    fig.savefig(fig_name)
    log.info(f"    saved figure '{fig_name}'")
    
    return (traj, )


@validators.validate_trajectory
def plot_rmsd_all_chains_one_subplot(
        traj,
        *args,
        chains="all",
        fig_name="plot_rmsd_all_chains_one_subplot.pdf",
        colors="None",
        **kwargs,
        ):
    """
    Plots the RMSDs of the selected chains in a single subplot.
    
    Parameters:
    
        - chains (str, opt, def: "all"): the chains to consider for
            plotting. Accepts integer slicer: 1:, :50, 2:30:2, etc...
        
        - colors (str of comma separated values): "blue,red,green".
        
        - fig_name (str, opt, def: "RMSDS_plot_combined.pdf"): the
            name of the output figure.
    """
    
    log.info("* Plotting chains RMSDs single subplot...")
    
    validators.validate_args(
        (
            (chains, str),
            (fig_name, str),
            (colors, str),
            ),
        func_name=plot_chains_single_subplot.__name__,
        )
    
    chain_list = pltcommons._get_chain_list(traj, chains)
    colors = pltcommons._get_colors(colors)
        
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    plt.tight_layout(rect=[0.05, 0.02, 0.995, 0.985])
    fig.suptitle("Chains' RMSDs", x=0.5, y=0.990, va="top", ha="center")
    
    x_range = np.arange(0, traj.n_frames)
    max_rmsds = 0
    for chain in chain_list:
        
        selector = f"chainid {chain}"
        slicer = traj.topology.select(selector)
        log.debug(
            f"performing on: selector {selector} "
            f"with slicer len {len(slicer)}"
            )
        
        if slicer.size == 0:
            log.debug("This chains is NOT found in trajectory... ignoring...")
            continue
        
        else:
            chain_traj = traj.atom_slice(slicer, inplace=False)
            _, rmsds = fromtraj.calc_rmsds(chain_traj)
        
        ax.plot(
            x_range,
            rmsds,
            label=chain,
            color=next(colors),
            alpha=0.7,
            )
        
        if rmsds.max() > max_rmsds:
            max_rmsds = rmsds.max()
    
    ax.set_xlim(x_range[0], x_range[-1])
    ax.set_ylim(0, max_rmsds + max_rmsds * 0.1)
    
    ax.set_xlabel("Frame Number", weight='bold')
    ax.set_ylabel("RMSDs", weight='bold')
    
    ax.grid(color='lightgrey', linestyle='-', linewidth=1, alpha=0.7)
    ax.legend(fontsize=6, loc=4)
    
    fig.savefig(fig_name)
    log.info(f"    saved figure '{fig_name}'")
    
    return (traj, )
