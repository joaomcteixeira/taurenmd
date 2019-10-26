"""
Calculate parameters.
"""

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
    
    log.info(S('for selection: {}'.format(selection)))
    log.info(S('for {} frames'.format(frame_slice)))
    
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
