"""
Functions that wrap around `MDAnalysis library`_.

Functions contained in this module operate with MDAnalysis (MDA)
functionalities, either by using MDA to access Molecular Dynamics
data or by receiving MDA data structures and parsing them in some way.

When using functions contained in this library you should `cite both`_
taurenmd and MDAnalysis.

.. _MDAnalysis library: https://www.mdanalysis.org
.. _cite both: https://taurenmd.readthedocs.io/en/latest/citing.html
"""
import MDAnalysis as mda

import taurenmd.core as tcore
from taurenmd import Path, log
from taurenmd.libs import libcli, libio
from taurenmd.logger import S, T


@libcli.add_reference(tcore.ref_mda)
def load_universe(topology, *trajectories):
    """
    Load MDAnalysis universe.

    Accepts MDAnalysis compatible `topology formats`_ and
    `trajectory formats`_. Read further on the `MDAnalysis Universe`_.

    .. _topology formats: https://www.mdanalysis.org/docs/documentation_pages/topology/init.html
    .. _trajectory formats: https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#supported-coordinate-formats
    .. _MDAnalysis Universe: https://www.mdanalysis.org/docs/documentation_pages/core/universe.html?highlight=universe#core-object-universe-mdanalysis-core-universe

    Examples
    --------
        
        >>> libmda.load_universe('topology.pdb', 'trajectory.dcd')

        >>> libmda.load_universe(
                'topology.pdb',
                'traj_part_1.xtc',
                'traj_part_2.xtc',
                Path('my', 'md', 'folder', 'traj_part_3.xtc'),
                )

    Parameters
    ----------
    topology : str or Path object
        Path to topology file.

    trajectories* : str of Path objects
        Paths to trajectory file(s). Trajectory files will be used
        sequentially to create the Universe.
    
    Return
    ------
    MDAnalysis Universe
    """  # noqa: E501 D412
    libio.report_input(topology, trajectories)
    universe = mda.Universe(
        Path(topology).str(),
        [Path(i).str() for i in trajectories],
        )
    report(universe)
    return universe


@libcli.add_reference(tcore.ref_mda)
def report(universe):
    """
    Report information about the Universe.
    
    Example
    -------

        >>> u = libmda.load_universe('topology.pdb', 'trajectory.xtc')
        >>> libmda.report(u)
       
    Parameters
    ----------
    universe : MDAnalysis Universe
        `MDAnalysis universe <https://www.mdanalysis.org/docs/documentation_pages/core/universe.html?highlight=universe#core-object-universe-mdanalysis-core-universe>`_.

    Returns
    -------
    None

    """  # noqa: E501
    log.info(T('Reporting'))
    log.info(S('number of frames: {}', len(universe.trajectory)))
    log.info(S('number of atoms: {}', len(universe.atoms)))


@libcli.add_reference(tcore.ref_mda)
def draw_atom_label_from_atom_group(atom_group):
    """
    Translate MDAnalysis Atom Group to list of strings for each atom.

    Strings represent each atom by SEGID.RESNUM|RESNAME.NAME,
    for example carbon alpha of Cys 18 of chain A would:
        
        >>> A.18Cys.CA

    This function is used by taurenmd for data representation
    purporses.

    Parameters
    ----------
    atom_group : Atom Group obj
       `MDAnalysis Atom group <https://www.mdanalysis.org/docs/documentation_pages/core/groups.html?highlight=atom%20group#MDAnalysis.core.groups.AtomGroup>`_.

    Returns
    -------
    list of strings
        Containing the atom string representation for each atom in
        the Atom Group.
    """  # noqa: E501
    labels = []
    for atom in atom_group:
        s = '{}.{}{}.{}'.format(
            atom.segment.segid,
            atom.residue.resnum,
            atom.resname.title(),
            atom.name,
            )
        labels.append(s)

    return labels
