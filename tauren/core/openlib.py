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

import os
import sys
import json
from collections import namedtuple

import mdtraj as md
import simtk.openmm.app as app

from tauren import logger
from tauren import system

log = logger.get_log(__name__)

_traj_details = """
*** Loaded ***

trajectory: {}
topology: {}

* details:

n_frames: {}
n_residues: {}
n_atmos: {}

time_step: {} ps or {} ns
total_time: {} ps or {} ns
"""


def load_json_config(config_path):
    """
    Loads configuration JSON file into a collections.namedtuple()
    
    Parameters:
    
        - config_path (str): path to .json file
    
    Returns:
    
        - namedtuple() with 'typename' = 'config and 'field_names'
            the JSON config keys.
    """
    
    assert config_path.endswith('.json'), \
        "config file should be JSON (.json) type"
    
    assert os.path.exists(config_path), \
        "'{}' file does not exists".format(config_path)
    
    assert os.path.isfile(config_path), \
        "'{}' is not a file".format(config_path)
    
    with open(config_path, 'r') as conf:
        config = json.load(conf)
    
    a = namedtuple("config", config.keys())
    config_tuple = a(**config)
    
    log.debug("<config_tuple>: {}".format(config_tuple))
    
    return config_tuple


def load_traj(traj_file, topo_file):
    """
    Loads MD trajectory.
    
    Parameters:
    
        - traj_file (str): trajectory file name (path)
            Formats allowed: ".xtc", ".nc", ".trr", ".h5", ".pdb",
                ".binpos", ".dcd"
        
        - topo_file (str): topology file name (path)
            Formats allowed: ".pdb" and ".cif"
    
    Returns:
    
        - mdtraj.Trajectory object
        (http://mdtraj.org/1.9.0/api/generated/mdtraj.Trajectory.html)
    """
    
    log.info("loading trajectory...")
    
    if not(os.path.exists(topo_file)):
        
        log.info("topology file: '{}' does not exists\n".format(topo_file))
        sys.exit(1)
    
    elif not(topo_file.endswith(system.topology_types)):
        log.info("topology file not suited\n")
        log.info("should have extention: '{}'\n".format(system.topology_types))
        sys.exit(1)
    
    if not(os.path.exists(traj_file)):
        log.info("trajectory file: '{}' does not exists\n".format(traj_file))
        sys.exit(1)
    
    elif not(traj_file.endswith(system.trajectory_types)):
        log.info("trajectory file not suited\n")
        log.info("should end with '{}'\n".format(system.trajectory_types))
        sys.exit(1)
    
    if topo_file.endswith(".cif"):
        structure = app.PDBxFile(topo_file)
        log.info("loaded topology with openmm.app.PDBxFile")
        topology = md.Topology.from_openmm(structure.topology)
        log.info("loaded topology with md.Topology.from_openmm")
    
    elif topo_file.endswith('.pdb'):
        topology = topo_file
    
    traj = md.load(traj_file, top=topology)
    
    info = _traj_details.format(
        traj_file,
        topo_file,
        traj.n_frames,
        traj.n_residues,
        traj.n_atoms,
        traj.timestep,
        traj.timestep / 1000,
        traj.time[-1],
        traj.time[-1] / 1000
        )
    
    log.info(info)
    
    return traj
