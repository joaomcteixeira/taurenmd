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

from tauren import tcomm, ttrans, tplot

trajectory_types = (".xtc", ".nc", ".trr", ".h5", ".pdb", ".binpos", ".dcd")
topology_types = (".pdb", ".cif")

actions_dict = {
    "remove_solvent": ttrans.transform.remove_solvent,
    "reduce_equidistant": ttrans.transform.reduce_equidistant,
    "try_mdtraj_image_molecules": ttrans.transform.mdtraj_image_molecules,
    "frames2PDB": tcomm.export.frames2PDB,
    "save_traj": tcomm.export.save_traj,
    "plot_rmsd_combined": tplot.rmsds.plot_rmsd_combined,
    "plot_rmsd_chain_per_subplot": tplot.rmsds.plot_rmsd_chain_per_subplot,
    "plot_rmsd_all_chains_one_subplot": tplot.rmsds.plot_rmsd_all_chains_one_subplot,
    }
