"""
# Report on trajectory characteristics.

## Example

     taurenmd report topology.pdb trajectory.dcd

## References

"""
import argparse
import functools

import taurenmd.core as tcore
from taurenmd import _BANNER, log
from taurenmd.libs import libcli, libmda
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += tcore.ref_mda

_help = 'Reports on trajectory details.'
_name = 'report'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)


def _ap():
    return ap


def main(topology, trajectories, **kwargs):
    """Execute main client logic."""
    log.info(T('reporting'))
    libmda.load_universe(topology, *trajectories)
    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
