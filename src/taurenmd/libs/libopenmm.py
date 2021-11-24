"""
Functions that wrap around `OpenMM library`_.

Read our `citing documentation`_ to understand how to cite multiple
libraries.

.. _citing documentation: https://taurenmd.readthedocs.io/en/latest/citing.html
.. _OpenMM library: http://openmm.org/
"""
import sys


try:
    from openmm.app.pdbxfile import PDBxFile
    SIMTK = True
except ImportError:
    SIMTK = False


from taurenmd import Path
from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libcli
from taurenmd.logger import S, T


@libcli.add_reference(tcore.ref_openmm)
def attempt_to_load_top_from_simtk(topology):
    """
    Load topology from SIMTK.

    Parameters
    ----------
    topology : str or Path

    Returns
    -------
    topology from mdtraj.Topology.from_openmm`

    Raises
    ------
    Dependency error from :func:`_log_simtkimport_error`, program
    halts.
    """
    topp = Path(topology)

    if topp.suffix == '.cif' and SIMTK:
        mol = PDBxFile(topp.str())
        return mol

    elif topp.suffix == '.cif' and not SIMTK:
        _log_simtkimport_error()

    else:
        raise ValueError(f'`topology` suffix is not CIF: {topp.suffix!r}')


def _log_simtkimport_error():
    msg = (
        "To use .cif files as Topologies taurenmd requires OpenMM, "
        "which is currently not installed. "
        "Please visit our Installation instruction at "
        "https://taurenmd.readthedocs.io/"
        )
    log.error(T('Dependency Error'))
    log.error(S(msg))
    sys.exit(0)
