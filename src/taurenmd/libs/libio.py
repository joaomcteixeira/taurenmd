"""
Handles input and output.
"""
from taurenmd import Path, log
from taurenmd.logger import S, T  # noqa: F401


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
