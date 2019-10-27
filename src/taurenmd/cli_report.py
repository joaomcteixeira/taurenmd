"""
Report on trajectory characteristics.
"""

from taurenmd import log
from taurenmd.libs import libcli


# import simtk.openmm as mm
# import simtk.openmm.app as app
# import mdtraj
import MDAnalysis as mda


ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

ap.add_argument(
    'topology',
    help='Topology file.',
    type=str,
    )

ap.add_argument(
    'trajectory',
    help='The trajectory file',
    nargs='?',
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
    log.info('Starting...')
    
    trj = libio.mda_load_universe(topology, *list(trajectory))
    return


if __name__ == '__main__':
    maincli()
