
from taurenmd import log
from taurenmd.logger import S


def frame_list(len_traj, start=None, stop=None, step=None, flist=None,):
    """
    Create frame integer list from a length and slice parameters.

    Parameters
    ----------
    start : int or None, optional
        The start index for the slice object.
        Defaults to ``None``.
    
    stop : int or None, optional
        The stop index for the slice object.
        Defaults to ``None``.

    step: int or None, optional
        the step index for the slice object.
        Defaults to ``None``.

    flist : list-like, or comma-separated string, optional
        The list of specific frames.
        Defaults to ``None``.
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
    

def _frame_slice(start=None, stop=None, step=None, ftuple=None):

    if isinstance(ftuple, (list, tuple)) and len(ftuple) == 3:
        sliceObject = slice(*[int(i) for i in ftuple])
    
    else:
        sliceObject = slice(start, stop, step)
    
    log.info(S('slicing: {}', sliceObject))
    return sliceObject


def evaluate_to_slice(*, value=None, start=None, stop=None, end=None):
    """
    Evaluate to slice.
    
    If any of ``start``, ``stop`` or ``step`` is given returns
    ``slice(start, stop, step)``. Otherwise tries to evaluate ``value``
    to its representative slice form.

    Examples
    --------

        >>> evalute_to_slice(value='1,100,2')
        >>> slice(1, 100, 2)
        
        >>> evaluate_to_slice(start=10)
        >>> slice(10, None, None)

        >>> evaluate_to_slice(value=(0, 50, 3))
        >>> slice(0, 50, 3)

        >>> evaluate_to_slice(value=(None, 100, None))
        >>> slice(None, 100, None)

        >>> #ATTENTION
        >>> evaluate_to_slice(value='10')
        >>> slice(10, None, None)
        >>> # this is different from slice(10)
        >>> #though
        >>> evaluate_to_slice(value=10)
        >>> slice(None, 10, None)

    Parameters
    ----------
    value : list, tuple, str, None or int
        A human readable value that can be parsed to a slice object
        intuitively.
        Defaults to ``None``.

    start : None or int
        The starting index of the slice (inclusive).
        Defaults to ``None``.

    stop : None or int
        The final index of the slice (exclusive).
        Defaults to ``None``.

    step : None or int
        Slice periodicity.
        Defaults to ``None``.
    
    Returns
    -------
    slice
        Python `slice object <https://docs.python.org/3/library/functions.html#slice>`_.

    Raises
    ------
    ValueError
        If slice can not be computed.
    """
    if any((start, stop, end)):
        return slice(start, stop, end)

    elif isinstance(value, (list, tuple)) and len(value) == 3:

        try:
            start = int(value[0])
        except (ValueError, TypeError):
            start = None

        try:
            stop = int(value[1])
        except (ValueError, TypeError):
            stop = None

        try:
            step = int(value[2])
        except (ValueError, TypeError):
            step = None

        return slice(start, stop, step)

    elif value is None:
        return slice(None, None, None)

    elif isinstance(value, str):
        if value.find(':') > -1:
            values = value.split(':')
        elif value.find(',') > -1:
            values = value.split(',')
        else:
            values = value.split()

        start, stop, step = [i if i else None for i in values]
        return slice(start, stop, step)
    
    elif isinstance(value, int):
        return slice(value)

    else:
        raise ValueError('slice object could not be generated')

