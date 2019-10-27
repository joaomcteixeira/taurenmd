"""
Input/Output functions using MDAnalysis.
"""
import MDAnalysis as mda

from taurenmd import log
from taurenmd.libs import libio
from taurenmd.logger import S, T


def mda_load_universe(topology, *trajectories):
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
    libio.report_input(topology, trajectories)

    universe = mda.Universe(topology, trajectories)

    mda_report(universe)
    return universe


def mda_report(universe):
    """Report information about the Universe."""
    log.info(T('Reporting'))
    log.info(S('number of frames: {}', len(universe.trajectory)))
    log.info(S('number of atoms: {}', len(universe.atoms)))
