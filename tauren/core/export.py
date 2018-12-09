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

from tauren import logger

log = logger.get_log(__name__)


def frames2PDB(traj, frames=None, suffix="_"):
    """
    Extracts trajectory frames to PDB files using <sufix> name.
    
    Parameters:
    
        - frames (opt, str): frame or range of frames to extract
            (INCLUSIVE), examples:
                "0", "0:", ":500", "10:50"
            Defaults to None -> extracts all frames
        
    """
    
    assert isinstance(suffix, str), "<suffix> must be str type."
    
    if frames:
        
        assert isinstance(frames, str), \
            "<suffix> must be str type."
        
        if frames.startswith(":"):
            try:
                frames_to_extract = list(range(0, int(frames[0:])))
            except ValueError:
                log.info("<frames> not valid: '{}'".format(frames))
                return
        
        elif frames.endswith(":"):
            try:
                frames_to_extrat = \
                    list(range(int(frames[:-1]), traj.n_frames + 1))
            except ValueError:
                log.info("<frames> not valid: '{}'".format(frames))
                return
        
        elif frames.find(","):
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
    
    else:
        frames_to_extract = range(traj.n_frames)
    
    if not(frames_to_extract):
        msg1 = "* ERROR * Couldn't identify the frames to extract: '{}'"
        log.info(msg1.format(frames_to_extract))
        
        msg2 = "* ERROR * Likely frames specified are outside the traj range?"
        log.info(msg2)

        return
    
    log.debug("<frames>: {}".format(frames))
    log.debug("<all_frames>: {}".format(all_frames))
    
    leading_zeros = len(str(traj.n_frames))
    pdb_name_fmt = suffix + "{:0>" + leading_zeros + "}"
    
    for frame in frames_to_extract:
        
        try:
            slice_ = traj.slice(frame, copy=True)
        
        except IndexError:
            msg = "* Frame '{}' does NOT exists in trajectory, ignoring...".\
                format(frame)
            log.info(msg)
            continue
        
        pdb_name = pdb_name_fmt.format(suffix, frame)
        slice_.save_pdb(pdb_name)
        log.info("* extracted {}".format(pdb_name))
    
    return
