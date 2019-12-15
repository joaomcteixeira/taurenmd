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


def draw_atom_label_from_atom_group(atom_group):
    """
    Translate `atom_group` to list of representing strings
    for each atom.
    """

    labels = []
    for atom in atom_group:
        s = '{}.{}{}.{}'.format(
            atom.segment.segid,
            atom.residue.resnum,
            atom.resname,
            atom.name,
            )
        labels.append(s)

    return labels
