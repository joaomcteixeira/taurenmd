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

from taurenmd.libs import libcli
import taurenmd.cli_traj2pdb as cli_traj2pdb
import taurenmd.cli_noSol as cli_noSol
import taurenmd.cli_imagemol as cli_imagemol


def load_args():
    
    ap = libcli.CustomParser(description=__doc__)

    subparsers = ap.add_subparsers(title='TAURENMD SUBROUTINES')

    ap_t2p = subparsers.add_parser(
        'traj2pdb',
        help='Converts traj to PDB.',
        parents=[cli_traj2pdb.ap],
        add_help=False,
        )
    ap_t2p.set_defaults(func=cli_traj2pdb.main)
    
    ap_noSol = subparsers.add_parser(
        'noSol',
        help='Removes solvent and extracts fist frame',
        parents=[cli_noSol.ap],
        add_help=False,
        )
    ap_noSol.set_defaults(func=cli_noSol.main)
    
    ap_imagemol = subparsers.add_parser(
        'imagemol',
        help='Attempts to image molecules.',
        parents=[cli_imagemol.ap],
        add_help=False,
        )
    ap_noSol.set_defaults(func=cli_imagemol.main)
    
    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    cmd = ap.parse_args()

    return cmd


def maincli():
    args = load_args()
    args.func(**vars(args))


if __name__ == '__main__':
    maincli()
