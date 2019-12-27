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
"""
import math

import numpy as np
from MDAnalysis.analysis.rms import RMSD as mdaRMSD
from MDAnalysis.analysis.rms import RMSF as mdaRMSF
from pyquaternion import Quaternion as Q

from taurenmd import log
from taurenmd.logger import S, T


def mda_rmsd_combined_chains(
        universe,
        frame_slice=None,
        selection='all',
        ref_frame=0,
        ):
    """
    Calculates RMSDs observed for a selection.
    
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

    frames : None or slice obj, optional
        The frames to consider. If ``None`` considers all frames.
        A slice object can be used for advanced slicing, for example:
        ``slice(None, 50, 2)``, slices from start to frame 50 (*exclusive*)
        every 2 frames.
        Defaults to `None`.

    selection : str, optional
        The selection upon which calculate the RMSDs.
        Defaults to `all`.

    ref_frames : int, optional
        The reference frame against which calculate the RMSDs.
        Defaults to 0.
    
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
        MDAnalysis related exceptions if something goes wrong with
        RMSD calculation.
    """
    log.info(T('Calculating RMSDs'))
    
    frame_slice = libutil.evaluate_to_slice(value=frame_slice)

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


def mda_rmsf(
        atom_group,
        frame_slice=None,
        ):
    """
    Calculates RMSFs.

    Uses `MDAnalysis RMSF <https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html?highlight=rmsd#MDAnalysis.analysis.rms.RMSF>`_.

    Parameters
    ----------
    atom_group : MDAnalysis Atom Group.
       `MDAnalysis Atom group <https://www.mdanalysis.org/docs/documentation_pages/core/groups.html?highlight=atom%20group#MDAnalysis.core.groups.AtomGroup>`_.

    frame_slice : any, optional
        Any argument that :py:func:`taurenmd.libs.libutil.evaluate_to_slice` can receive.
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
    """
    log.info(T('Calculating RMSFs'))
   
    frame_slice = libutil.evaluate_to_slice(value=frame_slice)
     
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
    Calculates the normal vector for the (p1, p2, p3) plane.

    Given 3 points, calculates the normal vector
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
    """
    cp = calc_plane_normal(p1, p2, p3)
    a, b, c = cp
    d = np.dot(cp, p3)
    return a, b, c, d


def calc_planes_angle(a1, b1, c1, a2, b2, c2):
    """
    Calculates the angle between two planes.

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

    Returns
    -------
    float
        The angle between plane 1 and plane 2 in degrees.
    """
    d = (a1 * a2 + b1 * b2 + c1 * c2)
    e1 = math.sqrt(a1 * a1 + b1 * b1 + c1 * c1)
    e2 = math.sqrt(a2 * a2 + b2 * b2 + c2 * c2)
    d = d / (e1 * e2)
    return math.degrees(math.acos(d))


def generate_quaternion_rotations(
        rotation_axis,
        rotating_vector,
        start=0,
        end=360,
        num=360 * 3,
        ):
    vec_u = Q(vector=rotating_vector).unit
    rot_u = Q(vector=rotation_axis).unit

    Q_rotated_tuples = []
    for angle in np.linspace(start, end, num=num):
        qq = Q(axis=rot_u.vector, degrees=angle)
        Q_rotated_tuples.append((qq, qq.rotate(vec_u)))

    return Q_rotated_tuples


def calc_minimum_Qdistances(rotation_tuples, reference_vector):
    
    ref_u = Q(vector=reference_vector).unit
    
    minimum = sorted(
        rotation_tuples,
        key=lambda x: Q.distance(x[1], ref_u),
        )
    return minimum[0][0]
