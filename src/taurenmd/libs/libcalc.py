"""
Calculate parameters.
"""
import math

import numpy as np
from MDAnalysis.analysis.rms import RMSD as mdaRMSD

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


def calc_plane_eq(p1, p2, p3):
    # http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
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
