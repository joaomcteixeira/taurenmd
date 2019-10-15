"""
Attempt image molecule with mdtraj.
"""
import argparse

import mdtraj

from taurenmd import log


def load_args():
    
    ap = argparse.ArgumentParser()

    ap.add_argument(
        'trajectory',
        help='The trajectory file.',
        )

    ap.add_argument(
        'topology',
        help='The topology file',
        )
    
    ap.add_argument(
        '-o',
        '--output',
        help='Output trajectory.',
        default='production_imaged.dcd',
        )
    
    cmd = ap.parse_args()
    return cmd


def main():
    cmd = load_args()
    main_script(**vars(cmd))


def main_script(trajectory, topology, output='production_imaged.dcd'):
    log.info('Starting...')
    trj = mdtraj.load(trajectory, top=topology)
    
    protein = trj.topology.guess_anchor_molecules()
    newtraj = trj.image_molecules(anchor_molecules=protein)
    
    # trj.image_molecules()
    newtraj.save(output)
    
    return


if __name__ == '__main__':
    main()
