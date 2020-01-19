"""
Functions that wrap around `MDTraj library`_.

Functions contained in this module operate with MDTraj functionalities,
wither by using MDTraj to access Molecular Dynamics data or by receiving
MDTraj data structures and parsing them in some way.

`Simtk OpenMM`_ is also used in some functions.

Read our `citing documentation`_ to understand how to cite multiple
libraries.

.. _MDTraj library: http://mdtraj.org/1.9.3/index.html
.. _citing documentation: https://taurenmd.readthedocs.io/en/latest/citing.html
.. _Simtk OpenMM: http://openmm.org/
"""
import sys

import mdtraj

import taurenmd.core as tcore
from taurenmd import Path, log
from taurenmd.libs import libcli, libio
from taurenmd.logger import S, T


try:
    import simtk.openmm.app as app
    SIMTK = True
except ImportError:
    SIMTK = False


def _log_simtkimport_error():
    msg = (
        "To use .cif files as Topologies taurenmd requires OpenMM, "
        "which is currently not installed. "
        "Please visit our Installation instruction at "
        "https://taurenmd.readthedocs.io/"
        )
    log.error(T('Dependency Error'))
    log.error(S(msg))
    sys.exit(0)


@libcli.add_reference(tcore.ref_mdt)
def load_traj(topology, trajectory):
    """
    Load trajectory with `MDTraj <http://mdtraj.org/1.9.3/index.html>`_.
    
    Uses `mdtraj.load <http://mdtraj.org/1.9.3/api/generated/mdtraj.load.html?highlight=load#mdtraj.load>`_.
    
    Example
    -------
        
        >>> libmdt.load_traj('bigtopology.cif', 'trajectory.dcd')

    Parameters
    ----------
    topology : str or Path
        Path to the topology file. Accepts MDTraj compatible `topology files <http://mdtraj.org/1.9.3/load_functions.html#trajectory-reference>`_. mmCIF format is loaded using `OpenMM <http://mdtraj.org/1.9.3/api/generated/mdtraj.Topology.html?highlight=from_openmm#mdtraj.Topology.from_openmm>`_.

    trajectory : str or Path
        Path to the trajectory file. Accepts MDTraj compatible `files <http://mdtraj.org/1.9.3/load_functions.html#trajectory-reference>`_

    Returns
    -------
    MDTraj trajectory
        `Trajectory object <http://mdtraj.org/1.9.3/api/generated/mdtraj.Trajectory.html#mdtraj-trajectory>`_.
    """  # noqa: E501
    libio.report_input(topology, trajectory)

    topp = Path(topology)
    if topp.suffix == '.cif' and SIMTK:
        mol = app.PDBxFile(topp.str())
        top = mdtraj.Topology.from_openmm(mol.topology)
    elif topp.suffix == '.cif' and not SIMTK:
        _log_simtkimport_error()
    else:
        top = topp.str()

    mdtrajectory = mdtraj.load(Path(trajectory).str(), top=top)

    return mdtrajectory


@libcli.add_reference(tcore.ref_mdt)
def imagemol_protocol1(traj):
    """Attempt to image molecules acting on the whole traj."""
    log.info(T('running reimage protocol #1'))
    log.info(S('finding molecules'))

    mols = traj.top.find_molecules()
    log.info(S('done'))
    
    log.info(T('reimaging'))
    reimaged = traj.image_molecules(
        inplace=False,
        anchor_molecules=mols[:1],
        other_molecules=mols[1:],
        )
    log.info(S('done'))
    return reimaged


@libcli.add_reference(tcore.ref_mdt)
def imagemol_protocol2(traj):
    """Attempt to image molecules frame by frame."""
    reimaged = []
    for frame in range(len(traj)):
        log.info(S('reimaging frame: {}', frame))
        
        mols = traj[frame].top.find_molecules()
    
        reimaged.append(
            traj[frame].image_molecules(
                inplace=False,
                anchor_molecules=mols[:1],
                other_molecules=mols[1:],
                )
            )

    log.info(S('concatenating traj frames'))
    # http://mdtraj.org/1.9.3/api/generated/mdtraj.join.html#mdtraj.join
    reimaged_traj = reimaged[0].join(reimaged[1:])

    return reimaged_traj
