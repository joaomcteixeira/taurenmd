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

from tauren import logger, system

log = logger.get_log(__name__)


def frames2PDB(traj, frames="all", suffix="_"):
    """
    Extracts trajectory frames to PDB files using <sufix> name.
    
    Parameters:
    
        - frames (opt, str): frame or range of frames to extract
            (INCLUSIVE), examples:
                "0", "0:", ":500", "10:50", "1,40,65,100"
            Defaults to "all" -> extracts all frames
        
    """
    
    assert isinstance(suffix, str), "<suffix> must be str type."
        
    assert isinstance(frames, str), \
        "<suffix> must be str type."
    
    if frames == "all":
        
        frames_to_extract = range(traj.n_frames)
    
    elif frames.startswith(":"):
        try:
            frames_to_extract = list(range(0, int(frames[1:]) + 1))
        except ValueError:
            log.info("<frames> not valid: '{}'".format(frames))
            return
    
    elif frames.endswith(":"):
        try:
            frames_to_extract = \
                list(range(int(frames[:-1]), traj.n_frames + 1))
        except ValueError:
            log.info("<frames> not valid: '{}'".format(frames))
            return
    
    elif bool(frames.find(",")):
        try:
            frames_to_extract = [int(f) for f in frames.split(',')]
        except ValueError:
            log.info("<frames> not valid: '{}'".format(frames))
            return
    else:
        try:
            frames_to_extract = list(int(frames))
        except ValueError:
            log.info("<frames> not valid: '{}'".format(frames))
            return
    
    if not(frames_to_extract):
        msg1 = "* ERROR * Couldn't identify the frames to extract: '{}'"
        log.info(msg1.format(frames_to_extract))
        
        msg2 = "* ERROR * Likely frames specified are outside the traj range?"
        log.info(msg2)

        return
    
    log.debug("<frames>: {}".format(frames))
    
    leading_zeros = str(len(str(traj.n_frames)))
    pdb_name_fmt = suffix + "{:0>" + leading_zeros + "}.pdb"
    print(pdb_name_fmt)
    
    for frame in frames_to_extract:
        
        try:
            slice_ = traj.slice(frame, copy=True)
        
        except IndexError:
            msg = "* Frame '{}' does NOT exists in trajectory, ignoring...".\
                format(frame)
            log.info(msg)
            continue
        
        pdb_name = pdb_name_fmt.format(frame)
        slice_.save_pdb(pdb_name)
        log.info("* extracted {}".format(pdb_name))
    
    return traj


def save_traj(traj, file_name="traj_output.dcd", overwrite=True):
    """
    Saves trajectory to <file_name>.
    Trajectory format is given by extension name.
    
    Uses mdtraj.Trajectory.save()
    
    Parameters:
    
        - file_name (str): name of the output trajectory file.
        
        - overwrite (bool): if file_name already exists, overwrites it.
    """
    log.info("* Exporting trajectory to: {}".format(file_name))
    
    if not(file_name.endswith(system.trajectory_types)):
        log.info("* ERROR * not a valid traj extension")
        log.info("* ERROR * should be '{}'".format(system.trajectory_types))
        log.info("* ignoring...")
        return
    
    traj.save(file_name, force_overwrite=overwrite)
    
    log.info("    ... saved")
    
    return traj
