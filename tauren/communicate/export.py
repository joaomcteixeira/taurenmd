"""
MODULE PROVIDES FUNCTIONS TO EXPORT TRAJECTORIES

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
import tauren.core as trncore

log = trncore.logger.get_log(__name__)


@trncore._decorators.validate_trajectory
def frames2PDB(
        traj,
        *args,
        frames="all",
        prefix="_",
        **kwargs,
        ):
    """
    Extracts trajectory frames to PDB files using <sufix> name.
    
    Parameters:
    
        - frames (opt, str): frame or range of frames to extract
            (INCLUSIVE), examples:
                "0", "0:", ":500", "10:50", "1,40,65,100"
            Defaults to "all" -> extracts all frames
        
        - prefix (opt, str): the prefix name for extracted PDBs.
            Defailts to "_"
        
    """
    
    if not isinstance(prefix, str):
        raise ValueError("<prefix> must be str type.")
    
    if not isinstance(frames, str):
        raise ValueError("<frames> must be str type.")
    
    log.debug(f"<frames>: {frames}")
    log.debug(f"<prefix>: {prefix}")
    
    frames_not_valid_str = f"<frames> input not valid: '{frames}'"
    frames_to_extract = 0
    
    possible_inputs = (
        "all" in frames,
        frames.startswith(":"),
        frames.endswith(":"),
        False if frames.find(",") < 0 else True,
        True,
        )
    
    dict_of_actions = {
        0: range(traj.n_frames),
        1: list(range(0, int(frames[1:]) + 1)),
        2: list(range(int(frames[:-1]), traj.n_frames + 1)),
        3: [int(f) for f in frames.split(',')],
        4: list(int(frames)),
        }
    
    which_action = possible_inputs.index(True)
    
    try:
        frames_to_extract = dict_of_actions[which_action]
    
    except ValueError as e:
        log.info(e)
        log.info(frames_not_valid_str)
        sys.exit(1)
    
    assert bool(frames_to_extract), \
        "* ERROR * Couldn't identify the frames to extract."
    
    leading_zeros = str(len(str(traj.n_frames)))
    pdb_name_fmt = prefix + "{:0>" + leading_zeros + "}.pdb"
    
    for frame in frames_to_extract:
        
        try:
            slice_ = traj.slice(frame, copy=True)
        
        except IndexError:
            log.info(
                f"* Frame '{frame}' does NOT exist in trajectory, ",
                "ignoring...",
                )
            continue
        
        pdb_name = pdb_name_fmt.format(frame)
        slice_.save_pdb(pdb_name)
        log.info(f"* extracted {pdb_name}")
    
    return (traj, )


@trncore._decorators.validate_trajectory
def save_traj(
        traj,
        *args,
        file_name="traj_output.dcd",
        overwrite=True,
        **kwargs,
        ):
    """
    Saves trajectory to <file_name>.
    Trajectory format is given by extension name.
    
    Uses mdtraj.Trajectory.save()
    
    Parameters:
    
        - file_name (str): name of the output trajectory file.
        
        - overwrite (bool): if file_name already exists, overwrites it.
    """
    log.info(f"* Exporting trajectory to: {file_name}")
    
    traj.save(file_name, force_overwrite=overwrite)
    
    log.info("    ... saved")
    
    return (traj, )
