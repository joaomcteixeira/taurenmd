"""
Client Report
=============

**Report on trajectory characteristics.**

Reports on various trajectory details.

**Example:**

    >>> taurenmd report topology.pdb trajectory.dcd

**References:**

"""
import argparse
import functools

from taurenmd import log
from taurenmd.libs import libcli, libmda
from taurenmd.logger import S, T

__doc__ += libcli.ref_mda

_help = 'Reports on trajectory details.'
_name = 'report'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

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
