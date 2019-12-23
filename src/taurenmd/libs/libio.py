"""
Handles input and output.
"""
from taurenmd import Path, log
from taurenmd.logger import S, T  # noqa: F401


def _get_number(path):
    return int(Path(path).stem.split('_')[-1])


def sort_numbered_input(*inputs):
    """
    Sort input paths to tail number.

    Input paths or strings should be formatted such:
        
        my_trajectory_#.dcd

    where # is the tag evaluated to int and input sorted accordingly.

    Returns
    -------
    list
        The sorted pathlist
    """
    return sorted(inputs, key=_get_number)
   

def report_input(topology, trajectory):
    log.info(S('loading trajectory: {}', trajectory))
    log.info(S('with topology: {}', topology))
    return


def get_stem(filename, ext=None):
    """
    Returns the stem name of a file name.
    
    Example:
        >>> get_stem('trajectory.dcd', ext='pdb')
        >>> trajectory.pdb
    
    Parameters
    ----------
    filename
        The file name.
    """
    
    output = Path(filename).stem

    if ext:
        output = Path(output).with_suffix('.' + ext.lstrip("."))

    return output


def mk_frame_path(input_path, frame=0, ext='.pdb'):
    """
    Create the path name for a frame.

    Given an input_path (usually the name of the trajectory), create
    the corresponding frame filename.

    Example:
        >>> mk_frame_path('traj_output.xtc')
        >>> traj_output_frame0.pdb
    """
    input_path = Path(input_path)
    top_output = Path(
        input_path.myparents(),
        '{}_frame{}'.format(input_path.stem, frame),
        ).with_suffix('.{}'.format(ext.lstrip('.')))

    return top_output


def save_to_file(xdata, ydata, fname, header='', fmt='{:.3}', delimiter=','):
    """
    Saves data to file.
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
