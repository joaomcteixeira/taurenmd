"""
Functions to calculate parameters from single trajectories.

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

import mdtraj as md

from tauren import tlog
from tauren._core import validators

log = tlog.get_log(__name__)


@validators.validate_trajectory
def calc_rmsds(traj, *args, ref_frame=0, **kwargs):
    """
    Calculates RMSDs from trajectory.
    """
    
    rmsds = md.rmsd(
        traj,
        traj,
        frame=ref_frame,
        parallel=True,
        precentered=False,
        )
    
    log.debug(
        f"<rmsds>: max {rmsds.max()},"
        f" min {rmsds.min()},"
        f" average {rmsds.mean()}"
        )
    
    return (traj, rmsds)
