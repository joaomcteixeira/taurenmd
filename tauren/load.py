"""
MODULE PROVIDES FUNCTIONS TO OPEN TRAJECTORIES

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
import functools
import json
from pathlib import Path
from collections import namedtuple

import mdtraj as md
import simtk.openmm.app as app

from tauren import logger
from tauren import tauren

log = logger.get_log(__name__)


def _validate_file_paths(func):
    """
    Validates paths. Paths should exists and be files.
    Raises Erros otherwise.
    
    Used in functions where all positional arguments are paths.
    """
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        
        for file_ in args:
            
            path_ = Path(file_)
            
            if not path_.exists():
                raise FileNotFoundError(f"'{path_}' does NOT exist.")
            
            if not path_.is_file():
                raise ValueError(f"'{path_}' is NOT a file.")
            
        return func(*args, **kwargs)
    
    return wrapper


@_validate_file_paths
def load_json_config(config_path):
    """
    Loads configuration JSON file into a collections.namedtuple()
    
    Parameters:
    
        - config_path (str): path to .json file
    
    Returns:
    
        - namedtuple() with 'typename' = 'config and 'field_names'
            the JSON config keys.
    """
    
    if not config_path.endswith('.json'):
        raise TypeError("config file should have '.json' extension.")
    
    with open(config_path, 'r') as conf:
        config = json.load(conf)
    
    # python dictionaries are OrderedDicts from +3.6
    a = namedtuple("config", config.keys())
    config_tuple = a(**config)
    
    log.debug(f"<config_tuple>: {config_tuple}")
    
    return config_tuple


@_validate_file_paths
def load_traj(
        traj_file,
        topo_file,
        traj_type="mdtraj",
        ):
    """
    Loads MD trajectory.
    
    Parameters
    ----------
    traj_file : str
        Trajectory file name (path)
            Formats allowed: ".xtc", ".nc", ".trr", ".h5", ".pdb",
                ".binpos", ".dcd"
        
    topo_file : str
        Topology file name (path)
            Formats allowed: ".pdb" and ".cif"
    
    traj_type : str
        Defaults to "mdtraj".
        The type of trajectory generated.
    
    Returns
    -------
    Tauren Trajectory
    """
    
    log.info("loading trajectory...")
    
    if topo_file.endswith(".cif"):
        structure = app.PDBxFile(topo_file)
        log.info("loaded topology with openmm.app.PDBxFile")
        topology = md.Topology.from_openmm(structure.topology)
        log.info("loaded topology with md.Topology.from_openmm")
    
    elif topo_file.endswith('.pdb'):
        topology = topo_file
    
    # Exceptions are handled directly by md.load()
    
    if traj_type == "mdtraj":
        
        traj = tauren.TaurenMDTraj(traj_file, topology)
    
    info = f"""
*** Loaded ***

trajectory: {traj_file}
topology: {topo_file}
"""
    
    log.info(info)
    
    traj.report()
    
    return traj
