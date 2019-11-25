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
import sys

import taurenmd.cli_angle as cli_angle
import taurenmd.cli_distances as cli_dist
import taurenmd.cli_fext as cli_fext
import taurenmd.cli_imagemol as cli_imagemol
import taurenmd.cli_noSol as cli_noSol
import taurenmd.cli_report as cli_report
import taurenmd.cli_rmsd as cli_rmsd
import taurenmd.cli_trajedit as cli_trajedit
from taurenmd import log
from taurenmd.libs import libcli


def load_args():
    
    ap = libcli.CustomParser(description=__doc__)

    subparsers = ap.add_subparsers(title='TAURENMD SUBROUTINES')

    ap_angle = subparsers.add_parser(
        'angle',
        help='Calculates the angle between a plane along the trajectory.',
        parents=[cli_angle.ap],
        add_help=False,
        )
    ap_angle.set_defaults(func=cli_angle.main)

    ap_dist = subparsers.add_parser(
        'dist',
        help='Calculates distances between geometric centres of selections',
        parents=[cli_dist.ap],
        add_help=False,
        )
    ap_dist.set_defaults(func=cli_dist.main)
    
    ap_imagemol = subparsers.add_parser(
        'imagemol',
        help='Attempts to image molecules.',
        parents=[cli_imagemol.ap],
        add_help=False,
        )
    ap_imagemol.set_defaults(func=cli_imagemol.main)
    
    ap_fext = subparsers.add_parser(
        'fext',
        help='Extracts single frames from trajectory.',
        parents=[cli_fext.ap],
        add_help=False,
        )
    ap_fext.set_defaults(func=cli_fext.main)
    
    ap_noSol = subparsers.add_parser(
        'noSol',
        help='Removes solvent and extracts fist frame',
        parents=[cli_noSol.ap],
        add_help=False,
        )
    ap_noSol.set_defaults(func=cli_noSol.main)
    
    ap_tedit = subparsers.add_parser(
        'trajedit',
        help='Converts traj to PDB.',
        parents=[cli_trajedit.ap],
        add_help=False,
        )
    ap_tedit.set_defaults(func=cli_trajedit.main)

    ap_report = subparsers.add_parser(
        'report',
        help='Reports on traj details.',
        parents=[cli_report.ap],
        add_help=False,
        )
    ap_report.set_defaults(func=cli_report.main)

    ap_rmsd = subparsers.add_parser(
        'rmsd',
        help='Calculates and plot RMSDs',
        parents=[cli_rmsd.ap],
        add_help=False,
        )
    ap_rmsd.set_defaults(func=cli_rmsd.main)

    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    cmd = ap.parse_args()

    return cmd


def maincli():
    args = load_args()
    log.debug(args)
    args.func(**vars(args))


if __name__ == '__main__':
    maincli()
