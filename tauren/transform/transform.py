"""
CONTAINS FUNCTIONS TO REDUCE TRAJECTORIES SIZE

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

import tauren.core as trncore

log = trncore.logger.get_log(__name__)


@trncore._decorators.validate_trajectory
def reduce_equidistant(
        traj,
        *args,
        step=1,
        **kwargs,
        ):
    """
    Reduces trajectory in equidistant frames separated by
    <step> frames.
    
    Parameters:
    
        - traj (MDTraj.Trajectory): the trajectory to modify
        
        - step (int): the number of steps between
        equidistant frames. Defaults to 1.
    
    Returns:
    
        - modified trajectory (MDTraj.Trajectory)
    """
    log.info("* Reducing Trajectory to equidistant frames...")
        
    if not isinstance(step, int):
        raise ValueError(f"<step> must be integer type: '{step}'")
    
    log.info(f"    received trajectory: {traj}")
    
    new_traj = traj[::step]
    
    log.info(f"    reduced trajectory: {new_traj}")
    
    return (new_traj, )


@trncore._decorators.validate_trajectory
def remove_solvent(
        traj,
        *args,
        exclude=None,
        **kwargs,
        ):
    """
    Removes solvent from Trajectory.
    
    Performs: MDTraj.Trajectory.remove_solvent()
    
    Parameters:
    
        - traj (MDTraj.Trajectory)
        
        - exclude: List of solvent residue names to retain
            in the new trajectory.
    
    Returns:
        
        - solventless trajectory (MDTraj.Trajectory)
    """
    log.info("* Removing solvent...")
    log.info(f"    received trajectory: {traj}")
    
    traj.remove_solvent(inplace=True, exclude=exclude)
    
    log.info(f"    solventless trajectory: {traj}")
    
    return (traj, )


@trncore._decorators.validate_trajectory
def mdtraj_image_molecules(
        traj,
        *args,
        anchor_molecules=None,
        other_molecules=None,
        sorted_bonds=None,
        make_whole=True,
        **kwargs,
        ):
    """
    Performs MDTraj.Trajectory.image_molecules, accepts same arguments.
    
    Returns:
    
        - modified trajectory
    """
    log.info("* Trying imaging molecules... this can take a while...")
    
    traj.image_molecules(
        inplace=True,
        anchor_molecules=anchor_molecules,
        other_molecules=other_molecules,
        sorted_bonds=sorted_bonds,
        make_whole=make_whole,
        )
    
    log.info("    completed.")
    
    return (traj, )
