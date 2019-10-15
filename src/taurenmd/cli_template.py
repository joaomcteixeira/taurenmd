"""
Does something.
"""
import argparse

from taurenmd import log


# import simtk.openmm as mm
# import simtk.openmm.app as app
# import mdtraj
# import MDAnalysis as mda


def load_args():
    """Load user arguments."""
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    
    ap.add_argument(
        'trajectory',
        help='The trajectory file',
        )
    
    ap.add_argument(
        'topology',
        help='Topology file.',
        type=str,
        )

    cmd = ap.parse_args()
    return cmd


def main():
    cmd = load_args()
    main_script(**vars(cmd))
    return


def main_script(trajectory):
    log.info('Starting...')
    return


if __name__ == '__main__':
    main()
