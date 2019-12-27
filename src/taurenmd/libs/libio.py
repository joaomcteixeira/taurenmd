"""Handle input and output general operations."""
from taurenmd import Path, log
from taurenmd.logger import S, T  # noqa: F401


def _get_number(path):
    return int(Path(path).stem.split('_')[-1])


def sort_numbered_input(*inputs):
    """
    Sort input paths to tail number.

    Input paths or strings should be formatted such::
        
        >>> my_trajectory_#.dcd

    Where ``#`` is the tag evaluated to ``int`` used to sort
    the input list.
    
    Parameters
    ----------
    *inputs : str of Paths
        Paths to files.

    Returns
    -------
    list
        The sorted pathlist.
    """
    return sorted(inputs, key=_get_number)
   

def report_input(topology, trajectory):
    """Report on topology and trajectory file paths used as input."""
    log.info(S('loading trajectory: {}', trajectory))
    log.info(S('with topology: {}', topology))
    return


def mk_frame_path(input_path, frame=0, ext='.pdb'):
    """
    Create the path name for a frame.

    Given an input_path (usually the name of the trajectory), create
    the corresponding frame file name.

    Example
    -------

        >>> mk_frame_path('traj_output.xtc')
        >>> traj_output_frame0.pdb

    Parameters
    ----------
    input_path : str or Path
        The file path. Normally, trajectory file path.

    frame : int, optional
        The frame to label the new path.

    ext : str, optional
        The returned file desired extension. Defaults to ``.pdb``.

    Returns
    -------
    Path
        The frame-labeled output path.
    """
    input_path = Path(input_path)
    top_output = Path(
        input_path.myparents(),
        '{}_frame{}'.format(input_path.stem, frame),
        ).with_suffix('.{}'.format(ext.lstrip('.')))

    return top_output


def export_data_to_file(
        xdata,
        ydata,
        fname,
        header='',
        fmt='{:.3}',
        delimiter=',',
        ):
    """
    Save data to file.

    Following the format:

        >>> # header
        >>> x1,ya1,yb1,yc1
        >>> x2,ya2,yb2,yc2
        >>> x3,ya3,yb3,yc3

    Where ``x`` are the elements in ``xdata``, and ``y*`` are the
    elements in the different ``ydata`` series.

    Parameters
    ----------
    xdata : list-like
        Contains the *x axis* data.

    ydata : list of list-like
        Contains the *y axis* data series.

    fname : str
        The output file path name.

    header : str
        The commented header of the file. Comment character, like ``#``
        is not placed, it should be already provided if desired.

    fmt : str
        The float format. Defaults to ``{:.3}``.

    delimiter : str
        The string delimiter between columns. Defaults to ``,``.
    """
    with open(fname, 'w') as fh:
        fh.write(header)
        for label, ydataseries in zip(
                xdata,
                zip(*ydata),
                ):

            fh.write('{}{}{}\n'.format(
                label,
                delimiter,
                delimiter.join(fmt.format(i) for i in ydataseries),
                ))
    return


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
    """  # noqa: E501
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
