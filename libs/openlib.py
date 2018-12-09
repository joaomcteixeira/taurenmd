"""
MODULE PROVIDES FUNCTIONS TO OPEN TRAJECTORIES
"""

import os
import sys

import logging
import numpy as np

import mdtraj as md
import simtk.openmm.app as app

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# create a file handler
debug_ = logging.FileHandler('md-tools.log')
debug_.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)

# create a logging format
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(filename)s:%(name)s:%(funcName)s:%(lineno)d - %(message)s')
debug_.setFormatter(formatter)

# add the handlers to the logger
log.addHandler(debug_)
log.addHandler(ch)

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
        
        sys.stderr.write("topology file: '{topo_file}' does not exists")
        sys.exit(1)
    
    if not(os.path.exists(traj_file)):
        sys.stderr.write("trajectory file: '{traj_file}' does not exists")
        sys.exit(1)
    
    if topo_file.endswith(".cif"):
        structure = app.PDBxFile(topo_file)
        log.info("loaded topology with openmm.app.PDBxFile")
        topology = md.Topology.from_openmm(structure.topology)
        log.info("loaded topology with md.Topology.from_openmm")
    
    elif topo_file.endswith('.pdf'):
        topology = topo_file
    
    else:
        sys.stderr.write("topology file not suited")
        sys.stderr.write("should have extention: 'pdb' or 'cif")
        sys.exit(1)
    
    trajectory_type = (".xtc", ".nc", ".trr", ".h5", ".pdb", ".binpos", ".dcd")
    
    if traj_file.endswith(trajectory_type):
        traj = md.load(traj_file, top=topology)
    
    else:
        sys.stderr.write("trajectory file not suited")
        sys.stderr.write("should end with {}".format(", ".join(trajectory_type)))
        sys.exit(1)
        
    info = \
"""
*** Loaded ***

trajectory: {}
topology: {}

* details:

n_frames: {}
n_residues: {}
n_atmos: {}

time_step: {} ps or {} ns
total_time: {} ps or {} ns
""".format(
        traj_file,
        topo_file,
        traj.n_frames,
        traj.n_residues,
        traj.n_atoms,
        traj.timestep,
        traj.timestep/1000,
        traj.time[-1],
        traj.time[-1]/1000
        )
    
    log.info(info)
    
    return traj
