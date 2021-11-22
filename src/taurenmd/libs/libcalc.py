"""
Calculate parameters from Molecular Dynamics data.

This module contains functions to calculate MD parameters such as:

#. RMSDs
#. RMSFs
#. plane angle variation
#. axes rotation decomposition

It contains also other functions that help on the calculation of the
desirables. Those functions are also avaible for independent use.

This library contains functions that operate on different Molecular
Dynamics data types. When special data types (MD analysis libraries)
are used, a prefix to the function name is used, and its docstring
explicitly refers to it.

When using these functions, you should always cite taurenmd together
with the other library(ies) used. `Read our citing reference page`_.

.. _Read our citing reference page: https://taurenmd.readthedocs.io/en/latest/citing.html
"""  # noqa: E501
import math

import numpy as np
from MDAnalysis.analysis.rms import RMSD as mdaRMSD
from MDAnalysis.analysis.rms import RMSF as mdaRMSF

from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libcli, libio
from taurenmd.logger import S, T


@libcli.add_reference(tcore.ref_mda)
def mda_rmsd(
        universe,
        frame_slice=None,
        selection='all',
        ref_frame=0,
        ):
    """
    Calculate RMSDs observed for a selection.

    Uses `MDAnalysis RMSD <https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html?highlight=rmsd#MDAnalysis.analysis.rms.RMSD>`_.

    Example
    -------

    * Calculate RMSDs observed for the whole system along the whole trajectory.

        >>> mda_rmsd(universe)

    * Calculate the RMSDs observed for selection `segid A` for every 10 frames.

        >>> mda_rmsd(universe, slice(0, None, 10), selection='segid A')

    Parameters
    ----------
    MDAnalysis Universe
        The MDAnalysis `universe <https://www.mdanalysis.org/docs/documentation_pages/core/universe.html?highlight=universe#core-object-universe-mdanalysis-core-universe>`_.

    frame_slice : any, optional
        The frames in the trajectory to consider.
        If ``None`` considers all frames.
        Accepts any argument that
        :py:func:`taurenmd.libs.libio.evaluate_to_slice` can receive.
        Defaults to ``None``.

    selection : str, optional
        The selection upon which calculate the RMSDs.
        Defaults to ``'all'``.

    ref_frames : int, optional
        The reference frame against which calculate the RMSDs.
        Defaults to ``0``.

    Returns
    -------
    Numpy Array
        The array containing the calculated RMSDs of shape (N,), where
        N is the number of frames.

    Raises
    ------
    ValueError
        If ``frame_slice`` is not ``None`` or ``slice`` object.

    MDAnalysis Exceptions
        Any exceptions that could come from MDAnalysis RMSF computation.
    """  # noqa: E501
    log.info(T('Calculating RMSDs'))

    frame_slice = libio.evaluate_to_slice(value=frame_slice)

    log.info(S('for selection: {}', selection))
    log.info(S('for trajectory slice: {}', frame_slice))

    R = mdaRMSD(
        universe,
        universe,
        select=selection,
        groupselection=None,
        ref_frame=ref_frame,
        verbose=False,
        )

    R.run(verbose=False)

    # rmsds[:, ii] = R.rmsd[:, 2][self._fslicer]

    return R.rmsd[frame_slice, 2]


@libcli.add_reference(tcore.ref_mda)
def mda_rmsf(
        atom_group,
        frame_slice=None,
        ):
    """
    Calculate RMSFs.

    Uses `MDAnalysis RMSF <https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html?highlight=rmsd#MDAnalysis.analysis.rms.RMSF>`_.

    Parameters
    ----------
    atom_group : MDAnalysis Atom Group.
       `MDAnalysis Atom group <https://www.mdanalysis.org/docs/documentation_pages/core/groups.html?highlight=atom%20group#MDAnalysis.core.groups.AtomGroup>`_.

    frame_slice : any, optional
        Any argument that :py:func:`taurenmd.libs.libio.evaluate_to_slice` can receive.
        Defaults to ``None``, considers all frames.

    Returns
    -------
    Numpy Array
        With the calculated RMSFs values, of shape (N,) where N are the
        frames sliced from ``frame_slice``.

    Raises
    ------
    MDAnalysis Exceptions
        Any exceptions that could come from MDAnalysis RMSF computation.
    """  # noqa: E501
    log.info(T('Calculating RMSFs'))

    frame_slice = libio.evaluate_to_slice(value=frame_slice)

    log.info(S('for trajectory slice: {}', frame_slice))

    R = mdaRMSF(
        atom_group,
        start=frame_slice.start,
        stop=frame_slice.stop,
        step=frame_slice.step,
        verbose=False,
        )

    R.run()

    return R.rmsf


def calc_plane_normal(p1, p2, p3):
    """
    Calculate the normal vector for the (p1, p2, p3) plane.

    Given 3 XYZ space coordinate points, calculates the normal vector
    of the plane defined by those points.

    Parameters
    ----------
    p1, p2, p3 : numpy.ndarray of shape (3,)
        The three 3D coordinate points that define the plane.

    Returns
    -------
    Numpy array of shape (3,)
        The normal vector to the (p1, p2, p3) plane. This vector
        is **NOT** an unitary vector.
    """
    v1 = p2 - p1
    v2 = p3 - p1
    return np.cross(v1, v2)  # normal to the plane


def calc_plane_eq(p1, p2, p3):
    """
    Calculate equation that defines the (p1, p2, p3) plane.

    `Further reading <http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/>`_.

    Parameters
    ----------
    p1, p2, p3 : numpy.ndarray of shape (3,)
        The three 3D coordinate points that define the plane.

    Returns
    -------
    tuple of length 4
        The four parameters (a, b, c, d) that defined the plane equation:

        .. math::

            ax + by + cz = d
    """  # noqa: E501
    cp = calc_plane_normal(p1, p2, p3)
    a, b, c = cp
    d = np.dot(cp, p3)
    return a, b, c, d


def calc_planes_angle(a1, b1, c1, a2, b2, c2, aunit='radians'):
    """
    Calculate the angle between two planes.

    Plane 1 is defined by a1, b1, c1 plane parameters,
    plane 2 is defined by a2, b2, c2, where:

    .. math::

        a1*x + b1*y + c1*z + d = 0

        a2*x + b2*y + c2*z + d = 0

    `Read further <from: https://www.geeksforgeeks.org/angle-between-two-planes-in-3d/>`_.

    Parameters
    ----------
    a1, b1, c1, a2, b2, c2 : float
        Plane parameters

    angle : str, optional
        ``degrees`` returns angle quantity in degrees, else returns
        in radians.

    Returns
    -------
    float
        The angle between plane 1 and plane 2.
    """  # noqa: E501
    d = (a1 * a2 + b1 * b2 + c1 * c2)
    e1 = math.sqrt(a1 * a1 + b1 * b1 + c1 * c1)
    e2 = math.sqrt(a2 * a2 + b2 * b2 + c2 * c2)
    d = d / (e1 * e2)
    angle = math.acos(d)
    if aunit == 'degrees':
        return math.degrees(angle)
    else:
        return angle


def calc_torsion_angles(coords):
    """
    Calculate torsion angles from sequential coordinates.

    Uses ``NumPy`` to compute angles in a vectorized fashion.
    Sign of the torsion angle is also calculated.

    Uses Prof. Azevedo implementation:
    https://azevedolab.net/resources/dihedral_angle.pdf

    Example
    -------
    Given the sequential coords that represent a dummy molecule of
    four atoms:

    >>> xyz = numpy.array([
    >>>     [0.06360, -0.79573, 1.21644],
    >>>     [-0.47370, -0.10913, 0.77737],
    >>>     [-1.75288, -0.51877, 1.33236],
    >>>     [-2.29018, 0.16783, 0.89329],
    >>>     ])

    A1---A2
           \
            \
            A3---A4

    Calculates the torsion angle in A2-A3 that would place A4 in respect
    to the plane (A1, A2, A3).

    Likewise, for a chain of N atoms A1, ..., An, calculates the torsion
    angles in (A2, A3) to (An-2, An-1). (A1, A2) and (An-1, An) do not
    have torsion angles.

    If coords represent a protein backbone consisting of N, CA, and C
    atoms and starting at the N-terminal, the torsion angles are given
    by the following slices to the resulting array:

    - phi (N-CA), [2::3]
    - psi (CA-N), [::3]
    - omega (N-C), [1::3]

    Parameters
    ----------
    coords : numpy.ndarray of shape (N>=4, 3)
        Where `N` is the number of atoms, must be equal or above 4.

    Returns
    -------
    numpy.ndarray of shape (N - 3,)
        The torsion angles in radians.
        If you want to convert those to degrees just apply
        ``np.degrees`` to the returned result.
    """
    # requires
    assert coords.shape[0] > 3
    assert coords.shape[1] == 3

    crds = coords.T

    # Yes, I always write explicit array indices! :-)
    q_vecs = crds[:, 1:] - crds[:, :-1]
    cross = np.cross(q_vecs[:, :-1], q_vecs[:, 1:], axis=0)
    unitary = cross / np.linalg.norm(cross, axis=0)

    # components
    # u0 comes handy to define because it fits u1
    u0 = unitary[:, :-1]

    # u1 is the unitary cross products of the second plane
    # that is the unitary q2xq3, obviously applied to the whole chain
    u1 = unitary[:, 1:]

    # u3 is the unitary of the bonds that have a torsion representation,
    # those are all but the first and the last
    u3 = q_vecs[:, 1:-1] / np.linalg.norm(q_vecs[:, 1:-1], axis=0)

    # u2
    # there is no need to further select dimensions for u2, those have
    # been already sliced in u1 and u3.
    u2 = np.cross(u3, u1, axis=0)

    # calculating cos and sin of the torsion angle
    # here we need to use the .T and np.diagonal trick to achieve
    # broadcasting along the whole coords chain
    # np.matmul is preferred to np.dot in this case
    # https://numpy.org/doc/stable/reference/generated/numpy.matmul.html
    cos_theta = np.diagonal(np.matmul(u0.T, u1))
    sin_theta = np.diagonal(np.matmul(u0.T, u2))

    # torsion angles
    return -np.arctan2(sin_theta, cos_theta)


def torsion_set(p1, p2, p3, p4_vecs):
    """
    Calculate torsion angles of set torwards three initial vectors.

    Same implementation as `calc_torsion_angles` but A1, A2, and A3 are
    fixed and A4 is an array of vectors. The torsion angles are
    calculated for A4s with respect to the fixed A1, A2, and A3.

    Example
    -------
    >>> torsion_set(
        [x1, y1, z1],
        [x2, y2, z2],
        [x3, y3, z3],
        [
            [x4_1, y4_1, z4_1],
            [x4_2, y4_2, z4_2],
            [x4_3, y4_3, z4_3],
            [x4_4, y4_4, z4_4],
            ...
            ],
        )

    Returns
    -------
    numpy.ndarray of shape (N,)
        The anges in radians. Where N is the length of A4.
    """
    assert p4_vecs.shape[1] == 3
    q1 = p2 - p1
    q2 = p3 - p2
    q3s = p4_vecs - p3
    assert q3s.shape == p4_vecs.shape

    q1x2 = np.cross(q1, q2)
    assert q1x2.shape == (3,)
    q2x3s = np.cross(q2, q3s)
    assert q2x3s.shape == q3s.shape

    n1 = q1x2 / np.linalg.norm(q1x2)
    assert n1.shape == (3,)
    n2 = q2x3s / np.linalg.norm(q2x3s)
    assert n2.shape == p4_vecs.shape, f'{n2.shape}, {p4_vecs.shape}'

    u1 = n2
    u3 = q2 / np.linalg.norm(q2)
    u2 = np.cross(u3, u1)
    assert u2.shape == p4_vecs.shape

    cos_theta = np.dot(n1, u1.T)
    sen_theta = np.dot(n1, u2.T)

    theta = -np.arctan2(sen_theta, cos_theta)
    assert theta.size == p4_vecs.shape[0]
    return theta
