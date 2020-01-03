"""Handle input and output general operations."""
import re

from taurenmd import Path, log
from taurenmd.logger import S, T  # noqa: F401


def get_number(path):
    """
    Extract tail number from path.
    
    Examples
    --------

        traj_1.dcd' -> 1
        traj_3.dcd' -> 3
        traj_1231.dcd' -> 1231
        traj_0011.dcd' -> 11
        traj_1_.dcd' -> 1
        traj_20200101_1.dcd -> 1

    Parameters
    ----------
    path : str or Path obj
        The path to evaluate.

    Returns
    -------
    int
        The tail integer of the path.
    """
    digit_re = re.compile('_\d+')
    stem = Path(path).stem
    number = re.findall('\d+', stem)[-1]
    return int(number)


def sort_numbered_input(*inputs):
    """
    Sort input paths to tail number.

    Sort criteria is provided by :py:func:`get_number`.

    Parameters
    ----------
    *inputs : str of Paths
        Paths to files.

    Returns
    -------
    list
        The sorted pathlist.
    """
    try:
        return sorted(inputs, key=get_number)
    except TypeError as err:
        log.exception(err)
        emsg = "Mind the packing *argument, input should be strings or Paths"
        raise ValueError(emsg)

   

def report_input(topology, trajectory):
    """Report on topology and trajectory file paths used as input."""
    log.info(S('loading trajectory: {}', trajectory))
    log.info(S('with topology: {}', topology))
    return


def add_prefix_to_path(ipath, prefix):
    """
    Add prefix to file path.
    
    Example
    -------

        >>> mk_frame_path('traj_output.xtc', prefix='my_prefix')
        >>> my_prefixtraj_output.xtc

    Mind the ``_`` is not placed automatically.

        >>> mk_frame_path('traj_output.xtc', prefix='my_prefix_')
        >>> my_prefix_traj_output.xtc

    Parameters
    ----------
    ipath : str or Path
        The file path to alter.

    prefix : str
        The complete prefix for the file name.

    Returns
    -------
    :py:func:`taurenmd.core.Path`
        The new file path.
    """
    ipath_ = _get_ipath(ipath)
    return Path(
        ipath_.myparents(),
        '{}{}'.format(prefix, ipath_.name),
        )


def add_suffix_to_path(ipath, suffix):
    """
    Add suffix to file path.
    
    If suffix has extention, updates the path extension, otherwise
    keeps the original extension.
    
    Examples
    --------

        >>> mk_frame_path('traj_output.xtc', suffix='my_suffix')
        >>> traj_outputmy_suffix.xtc

    Mind the underscore is not placed automatically:

        >>> mk_frame_path('traj_output.xtc', suffix='_my_suffix')
        >>> traj_output_my_suffix.xtc
   
    Updating extensions:

        >>> mk_frame_path('traj_output.xtc', suffix='_my_suffix.pdb')
        >>> traj_output_my_suffix.pdb

    Parameters
    ----------
    ipath : str or Path
        The file path to alter.

    suffix : str
        The complete suffix for the file name, extension should be
        included in the suffix, extension of the ``ipath`` is
        ignored.

    Returns
    -------
    :py:func:`taurenmd.core.Path`
        The new file path.
    """
    ipath_ = _get_ipath(ipath)
    result = Path(
        ipath_.myparents(),
        '{}{}'.format(ipath_.stem, suffix),
        )
    
    extension = Path(suffix).suffix or ipath_.suffix
    return result.with_suffix(extension)


def mk_frame_path(ipath, frame=0, ext='.pdb', leading=0, suffix=None):
    """
    Create the path name for a frame.

    Given an input_path ``ipath`` (usually the name of the trajectory),
    create the corresponding frame file name.

    Example
    -------

        >>> mk_frame_path('traj_output.xtc')
        >>> traj_output_frame0.pdb
       
        >>> mk_frame_path('traj_output.xtc', frame=4, leading=4)
        >>> traj_output_frame0004.pdb
   
    Parameters
    ----------
    ipath : str or Path
        The file path. Normally, trajectory file path.

    frame : int, optional
        The frame to label the new path.
        Defaults to ``0``.

    ext : str, optional
        The returned file desired extension.
        Defaults to ``.pdb``.

    leading : int
        The leading zeros to left append to the frame number.
        Defaults to ``0``.
    
    suffix : str
        Complete specifications of the desired suffix.
        If ``suffix`` is given, ``frame`` and ``ext`` and ``leading``
        are ignored and :py:func:`add_suffix_to_path` is used directly.

    Returns
    -------
    :py:func:`taurenmd.core.Path`
        The new file path.
    """
    if suffix:
        return add_suffix_to_path(ipath, suffix)
    else:
        ipath_ = _get_ipath(ipath)
        return Path(
            ipath_.myparents(),
            '{}_frame{}'.format(ipath_.stem, str(frame).zfill(leading)),
            ).with_suffix('.{}'.format(ext.lstrip('.')))


def _get_ipath(ipath):
    try:
        return Path(ipath)
    except TypeError as err:
        log.exception(err)
        return Path()


def parse_top_output(top_output, traj_output=None):
    """
    Parse different output definitions for topology output file name.
   
    If ``top_output`` startswith ``_`` uses :py:func:add_suffix_to_path.
    If ``top_output`` endswith ``_`` uses :py:func:add_prefix_to_path.
    Else: Return Path object of ``top_output``.


    Parameters
    ----------
    top_output : str or Path
        The string to evaluate.

    traj_output : str or Path
        The trajectory output file name. Return value depends on
        this parameters.

    Returns
    -------
    :py:func:`taurenmd.core.Path`
        The new topology file path.
    """
    if str(top_output).startswith('_'):
        return add_suffix_to_path(traj_output, suffix=top_output)

    elif str(top_output).endswith('_'):
        return add_prefix_to_path(traj_output, prefix=top_output)

    else:
        return Path(top_output)


def export_data_to_file(
        xdata,
        *ydata,
        fname='results.csv',
        header='',
        fmt='{:.3f}',
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
        fh.write('\n')
        for label, ydataseries in zip(
                xdata,
                zip(*ydata),
                ):

            fh.write('{}{}{}\n'.format(
                label,
                delimiter,
                delimiter.join(fmt.format(float(i)) for i in ydataseries),
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
    

def frame_slice(start=None, stop=None, step=None,):
    log.info(S('slicing: {}:{}:{}', start, stop, step))
    return slice(start, stop, step)


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
