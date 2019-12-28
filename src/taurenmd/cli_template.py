"""
Does something.
"""
import argparse

from taurenmd import CMDFILE, log
from taurenmd.libs import libcli, libio  # noqa: F401
from taurenmd.logger import S, T


ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
#libcli.add_trajectory_arg(ap)


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    libcli.save_command(CMDFILE, *sys.argv)
    main(**vars(cmd))
    return


def main(topology, trajectories, **kwargs):
#def main(topology, trajectory, **kwargs):
    log.info(T('starting'))
    
    #

    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()
