"""
Handles input and output.
"""

import MDAnalysis as mda

from taurenmd import log, Path
from taurenmd.logger import T, S


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
    universe = mda.Universe(top, traj)
    return universe


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
