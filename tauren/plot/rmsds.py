"""
CONTAINS GENERAL FUNCTIONS FOR PLOTTING

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

import numpy as np
import mdtraj as md
from matplotlib import pyplot as plt

import tauren.core as trncore

log = trncore.logger.get_log(__name__)


def _get_rmsds(traj):
    """
    Gets RMSD from traj.
    """
    
    rmsds = md.rmsd(traj, traj, frame=0, parallel=True, precentered=False)
    
    log.debug("<rmsds>: {}".format(rmsds))
    
    return rmsds


@trncore._decorators.validate_trajectory
def plot_overall_rmsd(
        traj,
        *args,
        color='blue',
        fig_name='overall_rmsds.pdf',
        **kwargs,
        ):
    """
    Plots a single plot with the combined RMSD for all chains.
    """
    
    log.info("* Plotting overall RMSDs...")
    
    rmsds = _get_rmsds(traj)
    
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    x_range = np.arange(1, traj.n_frames + 1)
    
    ax.plot(
        x_range,
        rmsds,
        label='all chains',
        color=color,
        alpha=0.7
        )
    
    ax.set_title('Overall RMSDs', weight='bold')
    ax.set_xlabel("Number of frames", weight='bold')
    ax.set_ylabel("RMSD", weight='bold')
    ax.grid(color='lightgrey', linestyle='-', linewidth=1, alpha=0.7)
    
    ax.set_xlim(x_range[0], x_range[-1])
    ax.set_ylim(0)
    
    ax.legend(fontsize=6, loc=4)
    
    fig.savefig(fig_name)
    log.info(f"    saved figure '{fig_name}'")
    
    return (traj, )
