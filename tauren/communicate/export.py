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

from tauren import tlog, tcommons
from tauren._core import validators

log = tlog.get_log(__name__)


@validators.validate_trajectory
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
    
    if frames == "all":
        frames_to_extract = range(traj.n_frames)
    else:
        frames_to_extract = tcommons.int_slicer(frames, traj.n_frames)
    
    leading_zeros = str(len(str(traj.n_frames)))
    pdb_name_fmt = prefix + "{:0>" + leading_zeros + "}.pdb"
    
    for frame in frames_to_extract:
        
        try:
            slice_ = traj.slice(frame, copy=True)
        
        except IndexError:
            log.exception(
                f"* Frame '{frame}' does NOT exist in trajectory, ",
                "ignoring...",
                )
            continue
        
        pdb_name = pdb_name_fmt.format(frame)
        slice_.save_pdb(pdb_name)
        log.info(f"* extracted {pdb_name}")
    
    return (traj, )


@validators.validate_trajectory
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
