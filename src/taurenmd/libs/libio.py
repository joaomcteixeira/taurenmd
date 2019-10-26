"""
Handles input and output.
"""

import MDAnalysis as mda
import mdtraj
import simtk.openmm.app as app

from taurenmd import Path, log
from taurenmd.logger import S, T  # noqa: F401


def mda_load_universe(top, traj):
    """
    Load MDAnalysis universe.

    Parameters
    ----------
    top
        Topology file.

    traj
        Trajectory file.
    
    Return
    ------
    MDAnalysis Universe
    """
    log.info(S(f'loading traj: {traj}'))
    log.info(S(f'loading top: {top}'))
    universe = mda.Universe(top, traj)
    log.info(S(S(f'number of atoms: {len(universe.atoms)}')))
    return universe


def mdtraj_load_traj(topology, traj):
    
    topp = Path(topology)
    if topp.suffix == '.cif':
        mol = app.PDBxFile(topp.str())
        top = mdtraj.Topology.from_openmm(mol.topology)
    else:
        top = topp.str()

    mdtrajectory = mdtraj.load(traj, top=top)

    return mdtrajectory


def get_stem(filename, ext=None):
    """
    Returns the stem name of a file name.
    
    Example:
        >>> get_stem('trajectory.dcd', ext='pdb')
        >>> trajectory.pdb
    
    Parameters
    ----------
    filename
        The file name.
    """
    
    output = Path(filename).stem

    if ext:
        output = Path(output).with_suffix(f'.{ext.lstrip(".")}')

    return output
