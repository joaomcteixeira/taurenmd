"""
CONTAINS GENERAL AND SYSTEM VARIABLES

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

from tauren.core import transform, export
from tauren.plot import general

trajectory_types = (".xtc", ".nc", ".trr", ".h5", ".pdb", ".binpos", ".dcd")
topology_types = (".pdb", ".cif")

actions_dict = {
    "remove_solvent": transform.remove_solvent,
    "reduce_equidistant": transform.reduce_equidistant,
    "try_mdtraj_image_molecules": transform.try_mdtraj_image_molecules,
    "frames2PDB": export.frames2PDB,
    "save_traj": export.save_traj,
    "plot_overall_rmsd": general.plot_overall_rmsd
    }
