"""
v0.1

Removes solvent from Trajectory

usage: python md_remove_solvent.py <TOPOLOGY> <TRAJECTORY INPUT> <TRAJECTORY OUTPUT>
        
        -TRAJECTORY, a MDTraj readable MD trajectory.
            According to: http://mdtraj.org/development/api/generated/mdtraj.Trajectory.html
        
        -TOPOLOGY, a MDTraj readable topology.
            According to: http://mdtraj.org/development/api/generated/mdtraj.Topology.html
    
examples:
    python md_frameextract.py file.pdb trajectory.dcd trajectory_out.dcd 

Author: {0} ({1})

"""

__author__ = "Joao M.C. Teixeira"
__email__ = "joaomcteixeira@gmail.com"

USAGE = __doc__.format(__author__, __email__)

import mdtraj as md
import sys
import os

def check_input(args):
    """ Checks input validity."""
    
    # confirms the number of arguments is correct
    if len(args) != 3:
        sys.stderr.write('* ERROR * args not passed correctly\n'.format(len(args)))
        sys.stderr.write(USAGE)
        sys.exit(1)
    
    trajectory_type = (".xtc", ".nc", ".trr", ".h5", ".pdb", ".binpos", ".dcd")
    for traj_file in args[1:]:
        if not traj_file.endswith(trajectory_type):
            sys.stderr.write('***ERROR: Not a valid trajectory file: ' + f + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)
    
    topology_type = (".pdb")
    if not args[0].endswith(topology_type):
        sys.stderr.write('Not a valid PDB file for MDTraj.load: ' + args[-1] + '\n')
        sys.stderr.write(USAGE)
        sys.exit(1)
    
    for f in args[:2]:
        if not os.path.isfile(f):
            sys.stderr.write('* ERROR * File not found: {}\n'.format(f))
            sys.stderr.write(USAGE)
            sys.exit(1)
    
    topology = args[0]
    traj_input = args[1]
    traj_output = args[2]
    
    return (topology, traj_input, traj_output)

def remove_solvent(traj):
    """Removes solvent"""
   
    print("Removing solvent...")
    traj.remove_solvent(inplace=True)
    print("* Solvent Removed")
    
    return traj

def save_trajectory(traj, traj_output_name):
    """
    Saves trajectory to file.
    """
    
    traj.save(traj_output_name)
    print("* Trajectory saved: {}".format(traj_output_name))
    
    return


if __name__ == "__main__":
    
    topology, traj_input, traj_output = check_input(sys.argv[1:])
    
    print("""
Options loaded:
- topology: {}
- traj_input: {}
- traj_output: {}
""".format(
        topology,
        traj_input,
        traj_output
        ))
    
    traj = md.load(traj_input, top=topology)
    print("Trajectory loaded: OK")
    
    print('How many frames has?     {}'.format(traj.n_frames))
    print('How many atoms has?      {}'.format(traj.n_atoms))
    print('How many residues has?   {}'.format(traj.n_residues))
    
    traj = remove_solvent(traj)
    
    save_trajectory(traj, traj_output)
    
