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
from MDAnalysis.analysis import align as mdaalign

from taurenmd import Path
from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libio, libopenmm, libutil
from taurenmd.logger import S, T


@tcore.add_reference(tcore.ref_mda)
def load_universe(
        topology,
        *trajectories,
        insort=False,
        **universe_args):
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

    insort : bool
        Whether to sort trajectory files by suffix number.
        See :func:`libio.sort_numbered_input`.

    universe_args : any
        Other arguments to be passed to `MDAnalysis Universe`.

    Return
    ------
    MDAnalysis Universe
    """  # noqa: E501 D412
    if insort:
        trajectories = libio.sort_numbered_input(*trajectories)

    libio.report_input(topology, trajectories)

    topo_path = Path(topology).str()
    traj_path = [Path(i).str() for i in trajectories],

    try:
        universe = mda.Universe(topo_path, traj_path, **universe_args)

    except ValueError:
        pdbx = libopenmm.attempt_to_load_top_from_simtk(topo_path)
        universe = mda.Universe(pdbx, traj_path, **universe_args)

    report(universe)
    return universe


@tcore.add_reference(tcore.ref_mda)
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
    segids = sorted(list(set(universe.atoms.segids)))
    info = {}
    for segid in segids:
        a = universe.select_atoms(f'segid {segid}')
        info[segid] = sorted(list(set(a.atoms.resids)))

    info_ = []
    for k, v in info.items():
        info_.append(str(S(
            f'segid {k} with {len(v)} residues '
            f'from {v[0]} to {v[-1]}'
            )))

    log.info(T('Reporting on universe'))
    log.info(S('number of frames: {}', len(universe.trajectory)))
    log.info(S('duration: {:.2f} ns', universe.trajectory[-1].time / 1000))

    ts2 = universe.trajectory[1].time
    ts1 = universe.trajectory[0].time
    log.info(S('timestep per frame: {:.2f} ns', (ts2 - ts1) / 1000))
    log.info(S('number of atoms: {}', len(universe.atoms)))
    log.info(S('components:\n{}', '\n'.join(info_)))


@tcore.add_reference(tcore.ref_mda)
def mdaalignto(universe, reference, selection='all'):
    """
    Align universe to reference.

    Uses `MDAnalysis.analysis.align.alignto <https://www.mdanalysis.org/docs/documentation_pages/analysis/align.html?highlight=alignto#MDAnalysis.analysis.align.alignto>`_.

    Parameters
    ----------
    universe, reference, selection
        Same as in ``MDAnalysis.analysis.align.alignto`` function.

    Raises
    ------
    ZeroDivisionError
        If selection gives empty selection.
    """  # noqa: E501
    try:
        mdaalign.alignto(universe, reference, select=selection)
    except (ValueError, ZeroDivisionError) as err:
        log.debug(err, exc_info=True)
        errmsg = (
            f'Could not perform alignment due to {err}, '
            'most likely the alignment selection does not match '
            'any possible selection in the system. You selection : '
            f'\'{selection}\'.'
            )
        log.info(errmsg)
        raise err


@tcore.add_reference(tcore.ref_mda)
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


@tcore.add_reference(tcore.ref_mda)
def convert_time_to_frame(x, dt, base_unit='ps'):
    """
    Convert a string `x` into a frame number based on given `dt`.

    If `x` does not contain any units its assumed to be a frame number
    already.

    Original function taken from `MDAnalysis.mdacli` project, from
    commit: https://github.com/MDAnalysis/mdacli/blob/15f6981df7b14ef5d52d64b56953d276291068ab/src/mdacli/utils.py#L25-L66
    The original function was modified internally without modifying its
    API.

    See also: https://github.com/MDAnalysis/mdacli/pull/81

    Parameters
    ----------
    x : str
        the input string
    dt : float
        the time step in ps

    Returns
    -------
    int
        frame number

    Raises
    ------
    ValueError
        The input does not contain any units but is not an integer.
    """  # noqa: E501
    # regex to split value and units while handling scientific input
    val, unit = libutil.split_time_unit(x)
    if unit != "":
        val = mda.units.convert(val, unit, base_unit)
        return int(val // dt)
    elif val % 1 != 0:  # the number is not int'able
        raise ValueError(
            "Only integers or time step combinations (`12ps`) "
            "are valid for frame selection."
            )
    else:
        return int(val)


def get_frame_list_from_slice(u, frame_slice):
    """
    Create a frame number list from a slice for a Universe.

    Parameters
    ----------
    u : mda.Universe
        The MDAnalysis universe.

    frame_slice : slice object

    Returns
    -------
    list of ints
        A list with the number of frames corresponding to that Universe
        and slice.
    """
    return list(range(len(u.trajectory))[frame_slice])


def create_x_data(u, xdata_in_time, frame_list):
    """."""
    if xdata_in_time:
        xdata = [
            mda.units.convert(
                u.trajectory[i].time - u.trajectory[0].time,
                'ps',
                xdata_in_time,
                )
            for i in frame_list
            ]
        xlabel = f'Time ({xdata_in_time})'
    else:
        xdata = frame_list
        xlabel = 'Frames'

    return xdata, xlabel
