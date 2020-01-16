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
from pyquaternion import Quaternion as Q

import taurenmd.core as tcore
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
    v1 = p3 - p1
    v2 = p2 - p1
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


@libcli.add_reference(tcore.ref_pyquaternion)
def generate_quaternion_rotations(
        rotation_axis,
        rotating_vector,
        start=0,
        end=360,
        num=360 * 3,
        endpoint=True,
        ):
    """
    Generate quaternion rotations of a vector around and axis.

    Rotates a `vector` around an `axis` for a series of angles.
    Rotation is performed using `Quaternion rotation <https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation>`_;
    namely, `pyquaterion.rotate <http://kieranwynn.github.io/pyquaternion/#rotation>`_.
    
    If you use this function you should cite
    `PyQuaternion package <http://kieranwynn.github.io/pyquaternion/>`_.

    Parameters
    ----------
    rotation_axis : tuple, list or array-like
        The XYZ coordinates of the rotation axis. This is the axis
        around which the vector will be rotated.

    rotation_vector : tuple, list or array-like
        The XYZ coordinates of the vector to rotate.

    start : int
        The starting rotation angle in degrees in the series.
        Defaults to ``0``.

    end : int
        The final rotation angle in degrees in the series.
        Defaults to ``360``.

    num : int
        The number of rotation steps in the angle rotation series.
        This is provided by ``np.linspace(start, end, num=num)``.

    Returns
    -------
    list of 2 unit tuples
        Where the list indexes follow the progression along
        ``np.linspace(start, end, num=num)`` and tuple elements
        are the rotation quatertion used to rotate ``rotation_vector``
        and the resulting rotated vector in its unitary form.
    """  # noqa: E501
    vec_u = Q(vector=rotating_vector).unit
    rot_u = Q(vector=rotation_axis).unit

    Q_rotated_tuples = []
    for angle in np.linspace(start, end, num=num, endpoint=endpoint):
        qq = Q(axis=rot_u.vector, degrees=angle)
        Q_rotated_tuples.append((qq, qq.rotate(vec_u)))

    return Q_rotated_tuples


@libcli.add_reference(tcore.ref_pyquaternion)
def sort_by_minimum_Qdistances(rotation_tuples, reference_vector):
    """
    Sort a list of quaternion rotation tuples.

    By its distance to a ``reference_vector``. Designed to receive
    the output of :py:func:`generate_quaternion_rotations`.

    If you use this function you should cite
    `PyQuaternion package <http://kieranwynn.github.io/pyquaternion/>`_.

    Parameters
    ----------
    rotation_tuples : list of tuples
        List of tuples containing Rotation quaternions and resulting
        rotated vector. In other words, the output from
        :py:func:`generate_quaternion_rotations`.

    reference_vector : tuple, list or array-like
        The XYZ 3D coordinates of the reference vector. The quaternion
        distance between vectors in ``rotation_tuples`` and this
        vector will be computed.
    
    Returns
    -------
    list
        The list sorted according to the creterion.
    """
    ref_u = Q(vector=reference_vector).unit
    
    minimum = sorted(
        rotation_tuples,
        key=lambda x: math.degrees(Q.distance(x[1], ref_u)) % 360,
        )
    return minimum
