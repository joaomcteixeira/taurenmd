"""
Does something.
"""
import argparse

from taurenmd import log
from taurenmd.libs import libcli, libio  # noqa: F401
from taurenmd.logger import S, T


ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    # formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    'topology',
    help='Topology file.',
    type=str,
    )

ap.add_argument(
    'trajectory',
    help=(
        'Trajectory files. If given, multiple trajectories will be'
        'contactenated by order.'
        ),
    nargs='+',
    )


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main(**vars(cmd))
    return


def main(topology, trajectory, **kwargs):
    log.info(T('starting'))
    
    #

    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()
