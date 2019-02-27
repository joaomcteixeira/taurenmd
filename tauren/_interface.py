"""
Contains variables to manage user-software interface.

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

from tauren import produce

actions_dict = {
    "remove_solvent": lambda x, y: x.remove_solvent(**y),
    "slice": lambda x, y: x.slice(**y),
    "try_image_molecules": lambda x, y: x.image_molecules(**y),
    "frames2file": lambda x, y: x.frames2file(**y),
    "save_traj": lambda x, y: x.save_traj(**y),
    "produce_rmsds_combined_chains":
        lambda x, y: produce.rmsds_combined_chains(x, **y),
    "produce_rmsds_separated_chains":
        lambda x, y: produce.rmsds_separated_chains(x, **y),
    }
