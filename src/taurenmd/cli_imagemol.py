"""
Attempt image molecule with mdtraj.
"""
import argparse

import mdtraj

from taurenmd import log
from taurenmd.libs import libcli


ap = libcli.CustomParser()

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


def load_args():
    cmd = ap.parse_args()
    return cmd


def main_cli():
    cmd = load_args()
    main_script(**vars(cmd))


def main(
        trajectory,
        topology,
        output='production_imaged.dcd',
        ):

    log.info('Starting...')
    
    trj = mdtraj.load(trajectory, top=topology)
   
    # use largest part as anchor
    mols = trj.top.find_molecules()
    reimaged = trj.image_molecules(
        inplace=False,
        anchor_molecules=mols[:1],
        other_molecules=mols[1:],
        )
    
    # trj.image_molecules()
    reimaged.save(output)
    
    return


if __name__ == '__main__':
    main_cli()
