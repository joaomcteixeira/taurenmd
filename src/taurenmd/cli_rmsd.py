"""
Does something.
"""
import argparse


# import simtk.openmm as mm
# import simtk.openmm.app as app
# import mdtraj
import MDAnalysis as mda

from taurenmd import log
from taurenmd.libs import libcli, libio, libcalc


_SELECTION = 'all'
_REF_FRAME = 0
_FRAMES = 'all'


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
    )

ap.add_argument(
    '-f',
    '--frames',
    help='The frames to consider.',
    nargs='+',
    default=_FRAMES,
    )

ap.add_argument(
    '-s',
    '--selection',
    help='The atom selection',
    type=str,
    default=_SELECTION,
    )

ap.add_argument(
    '-r',
    '--ref-frame',
    help='The refernece frame.',
    type=int,
    default=_REF_FRAME,
    )


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main_script(**vars(cmd))
    return


def main(
        topology,
        trajectory,
        frames=_FRAMES,
        ref_frame=_REF_FRAME,
        selection=_SELECTION,
        **kargs,
        ):
    log.info('Starting...')
    
    u = libio.mda_load_universe(topology, trajectory)

    rmsds_combined = libcalc.mda_rmsd_combined_chains(
        u,
        frames=frames,
        selection=selection,
        ref_frame=ref_frame,
        )


    return


if __name__ == '__main__':
    maincli()
