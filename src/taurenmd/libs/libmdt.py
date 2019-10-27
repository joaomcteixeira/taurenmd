"""
Input/Output functions for MDTraj.
"""
import mdtraj
import simtk.openmm.app as app

from taurenmd import Path
from taurenmd.libs import libio


def mdtraj_load_traj(topology, trajectory):
    """
    Loads trajectory with MDTraj.

    Parameters
    ----------
    topology
        The topology file.

    traj
        The trajectory file.

    Returns
    -------
    MDTraj trajectory.
    """
    libio.report_input(topology, trajectory)

    topp = Path(topology)
    if topp.suffix == '.cif':
        mol = app.PDBxFile(topp.str())
        top = mdtraj.Topology.from_openmm(mol.topology)
    else:
        top = topp.str()

    mdtrajectory = mdtraj.load(trajectory, top=top)

    return mdtrajectory
