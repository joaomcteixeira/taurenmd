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
