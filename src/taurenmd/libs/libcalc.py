"""
Calculate parameters.
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
    Calculates combined RMSDs.

    Combined RMSDs considers the selection as a whole.

    Parameters
    ----------
    universe
        The MDAnalysis universe.

    frames : str or tuple, optional
        The frames to consider.
        If 'all' considers all frames.
        Otherwise a tuple is needed of the format (start, end, step),
        to be used as a slice object.
        Defaults to `all`.

    selection : str, optional
        The selection upon which calculate the RMSDs.
        Defaults to `all`.

    ref_frames : int, optional
        Defaults to 0.
    
    Returns
    -------
    Numpy Array
        The array containing the calculated RMSDs.
    """
    log.info(T('Calculating RMSDs'))
    
    if frame_slice is None:
        frame_slice = slice(None, None, None)
    elif isinstance(frame_slice, slice):
        pass
    else:
        raise ValueError('frame_slice should be None or slice type')
    
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
    Calculates combined RMSDs.

    Combined RMSDs considers the selection as a whole.

    Parameters
    ----------
    universe
        The MDAnalysis universe.

    frames : str or tuple, optional
        The frames to consider.
        If 'all' considers all frames.
        Otherwise a tuple is needed of the format (start, end, step),
        to be used as a slice object.
        Defaults to `all`.

    selection : str, optional
        The selection upon which calculate the RMSDs.
        Defaults to `all`.

    ref_frames : int, optional
        Defaults to 0.
    
    Returns
    -------
    Numpy Array
        The array containing the calculated RMSDs.
    """
    log.info(T('Calculating RMSFs'))
    
    if frame_slice is None:
        frame_slice = slice(None, None, None)
    elif isinstance(frame_slice, slice):
        pass
    else:
        raise ValueError('frame_slice should be NoneType or Slice')
    
    log.info(S('for trajectory slice: {}', frame_slice))
    
    R = mdaRMSF(
        atom_group,
        start=frame_slice.start,
        stop=frame_slice.stop,
        step=frame_slice.step,
        verbose=False)
    
    R.run()
    
    return R.rmsf


def calc_plane_normal(p1, p2, p3):
    """
    Given 3 points, calculates the normal vector
    of the plane defined by those points.
    """
    v1 = p3 - p1
    v2 = p2 - p1
    return np.cross(v1, v2)  # normal to the plane


def calc_plane_eq(p1, p2, p3):
    # http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
    cp = calc_plane_normal(p1, p2, p3)
    a, b, c = cp
    d = np.dot(cp, p3)
    return a, b, c, d


def calc_angle(a1, b1, c1, a2, b2, c2):
    # from: https://www.geeksforgeeks.org/angle-between-two-planes-in-3d/
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
