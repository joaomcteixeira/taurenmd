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

import mdtraj as md

def equidistant_frames(traj, step):
    """
    Reduces trajectory in equidistant frames separated by
    <step> frames.
    
    Parameters:
    
        - traj (MDTraj.Trajectory): the trajectory to modify
        
        - step (int): the number of steps between
        equidistant frames.
    
    Returns:
    
        - modified trajectory (MDTraj.Trajectory)
    """
    
    assert isinstance(traj, md.Trajectory), "Not valid trajectory type"
    assert isinstance(step, int), "<step> must be integer type"
    
    return traj[::step]
