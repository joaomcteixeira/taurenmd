from taurenmd import log
from taurenmd.logger import S


def frame_list(len_traj, start=None, stop=None, step=None, flist=None,):
    """
    Createa frame list based on the length of trajectory.

    Parameters
    ----------
    start : int or None, optional
        The start index for the slice object.
    
    stop : int or None, optional
        The end index for the slice object.

    step: int or None, optional
        the step index for the slice object.

    flist : list-like, or comma-separated string, optional
        The list of specific frames.
        Defaults to None.
    """
    if any((start, stop, step)):
        return range(len_traj)[slice(start, stop, step)]

    elif flist:
        try:
            return [int(i) for i in flist.split(',')]
        except AttributeError:
            return [int(i) for i in flist]
        except Exception:
            raise ValueError(
                'Cannot generate list of frames from {}'.format(flist)
                ) from None
    else:
        return range(len_traj)
    

def frame_slice(start=None, stop=None, step=None, ftuple=None):

    if isinstance(ftuple, (list, tuple)) and len(ftuple) == 3:
        sliceObject = slice(*[int(i) for i in ftuple])
    
    else:
        sliceObject = slice(start, stop, step)
    
    log.info(S('slicing: {}', sliceObject))
    return sliceObject
