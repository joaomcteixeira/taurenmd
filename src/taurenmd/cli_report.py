"""
Report on trajectory characteristics.
"""
import argparse
import functools

from taurenmd import log
from taurenmd.libs import libcli, libmda
from taurenmd.logger import S, T

_help = 'Reports on trajectory details.'
_name = 'report'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)


def main(topology, trajectories, **kwargs):
    log.info(T('reporting'))
    libmda.load_universe(topology, *trajectories)
    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
