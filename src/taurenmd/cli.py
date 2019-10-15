"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mtaurenmd` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``taurenmd.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``taurenmd.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import argparse
import sys

from taurenmd.cli_traj2pdb import ap as ap_traj2pdb
from taurenmd.cli_traj2pdb import main_script as ms_traj2pdb


def load_args():
    
    ap = argparse.ArgumentParser(description=__doc__)

    subparsers = ap.add_subparsers(title='TAURENMD SUBROUTINES')

    ap_t2p = subparsers.add_parser(
        'traj2pdb',
        help='Converts traj to PDB.',
        parents=[ap_traj2pdb],
        add_help=False,
        )
    ap_t2p.set_defaults(func=ms_traj2pdb)

    if len(sys.argv) < 3:
        ap.print_help()
        ap.exit()

    cmd = ap.parse_args()

    return cmd


def main(args=None):
    args = load_args()
    args.func(**vars(args))
