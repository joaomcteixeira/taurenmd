"""
v0.1

Extracts frames from a MD trajectory using MDTraj.

usage: python md_frameextract.py -[NUM] -[OPT] <TRAJECTORY> <TOPOLOGY>

    OPTIONS:
        -NUM:
            -if integer, extract that frame.
            -if "all", extracts every frame.
        
        -OPT:
            -e: removes solvent molecules.
        
        -TRAJECTORY, a MDTraj readable MD trajectory.
            According to: http://mdtraj.org/development/api/generated/mdtraj.Trajectory.html
        
        -TOPOLOGY, a MDTraj readable topology.
            According to: http://mdtraj.org/development/api/generated/mdtraj.Topology.html
    
examples:
    python md_frameextract.py -1 trajectory.dcd file.pdb
    python md_frameextract.py -all trajectory.dcd file.pdb
    python md_frameextract.py -all -e trajectory.dcd file.pdb

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
    if not (3 <= len(args) <= 4):
        sys.stderr.write('*** ERROR: Too many options given: {}\n'.format(len(args)))
        sys.stderr.write(USAGE)
        sys.exit(1)
    
    if args[0].lstrip('-') == 'all':
        frames = 'all'
    else:
        try:
            frames = int(args[0].lstrip("-"))
            
        except ValueError:
            sys.stderr.write('***ERROR: Not a valid option: {}\n'.format(args[0]))
            sys.stderr.write(USAGE)
            sys.exit(1)
        
    if len(args) == 4:
        if args[1].lstrip('-') == 'e':
            options = ['e']
        
        else:
            sys.stderr.write('***ERROR: Not a valid option: {}\n'.format(args[1]))
            sys.stderr.write(USAGE)
            sys.exit(1)
    else:
        options = []
    
    for f in args[-2:]:
        if not os.path.isfile(f):
            sys.stderr.write('***ERROR: File not found: {}.\n'.format(f))
            sys.stderr.write(USAGE)
            sys.exit(1)
    
    trajectory_type = (".xtc", ".nc", ".trr", ".h5", ".pdb", ".binpos", ".dcd")
    if not args[-2].endswith(trajectory_type):
        sys.stderr.write('***ERROR: Not a valid trajectory file: ' + f + '\n')
        sys.stderr.write(USAGE)
        sys.exit(1)
    
    topology_type = (".pdb")
    if not args[-1].endswith(topology_type):
        sys.stderr.write('Not a valid PDB file for MDTraj.load: ' + args[-1] + '\n')
        sys.stderr.write(USAGE)
        sys.exit(1)
    
    return (frames, options, args[-2], args[-1])

def extract_frame(traj, i):
    """Extracts frame to PDB file."""
    
    try:
        slice_ = traj.slice(int(i), copy=True)
    
    except IndexError:
        msg = "*** You required to extract frame number {}. But this is not \
part of the trajectory. Please choose a frame wisely\n".format(i)
        sys.stderr.write(msg)
        sys.exit(1)
   
    slice_.save_pdb('_{:0>5}.pdb'.format(i))
    print("extracted _{:0>5}.pdb".format(i))
    
    return


if __name__ == "__main__":
    
    frames, options, trajf, topoh = check_input(sys.argv[1:])
    print("""
Options loaded:
- frames: {}
- options: {}
- Trajecotory file: {}
- Topology file: {}
""".format(
        frames,
        options,
        trajf,
        topoh
        ))
    
    traj = md.load(trajf, top=topoh)
    print("Trajectory loaded: OK")
    
    print('How many frames has?     {}'.format(traj.n_frames))
    print('How many atoms has?      {}'.format(traj.n_atoms))
    print('How many residues has?   {}'.format(traj.n_residues))
    
    if "e" in options:
        print("Removing solvent...")
        traj.remove_solvent(inplace=True)
        print("* Solvent Removed")
    
    print("Extracting frames...")
    if frames == 'all':
        for i in range(traj.n_frames):
            extract_frame(traj, i)
    else:
        extract_frame(traj, frames)
    
