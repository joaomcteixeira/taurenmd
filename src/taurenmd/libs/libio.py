"""
Handles input and output.
"""

import MDAnalysis as mda


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
