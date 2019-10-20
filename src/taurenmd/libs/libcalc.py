"""
Calculate parameters.
"""

import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD as mdaRMSD

from taurenmd import log
from taurenmd.logger import T, S


def mda_rmsd_combined_chains(
        universe,
        frames='all',
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
    
    if frames == 'all':
        frames = slice(None, None, None)
    else:
        frames = slice(*frames)
    
    log.info(S('for selection: {}'.format(selection)))
    log.info(S('for {} frames'.format(frames)))
    
    R = mdaRMSD(
        universe,
        universe,
        select=selection,
        groupselection=None,
        ref_frame=ref_frame,
        verbose=False,
        #superposition=False,
        #center=False,
        )
    
    R.run(verbose=False)
    
    # rmsds[:, ii] = R.rmsd[:, 2][self._fslicer]
    
    return R.rmsd[frames, 2]
