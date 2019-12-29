"""
Does something.
"""
import argparse
import functools

from taurenmd import log
from taurenmd.libs import libcli, libio  # noqa: F401
from taurenmd.logger import S, T


ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
#libcli.add_trajectory_arg(ap)


def main(topology, trajectories, **kwargs):
#def main(topology, trajectory, **kwargs):
    log.info(T('starting'))
    
    #

    log.info(S('done'))
    return

maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
