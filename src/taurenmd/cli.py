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
import taurenmd.cli_rmsf as cli_rmsf
import taurenmd.cli_rotations as cli_rot
import taurenmd.cli_trajedit as cli_trajedit
from taurenmd import CMDFILE, log
from taurenmd.libs import libcli


ap = libcli.CustomParser(description=__doc__)
subparsers = ap.add_subparsers(title='TAURENMD SUBROUTINES')

libcli.add_subparser(subparsers, cli_angle)
libcli.add_subparser(subparsers, cli_dist)
libcli.add_subparser(subparsers, cli_fext)
libcli.add_subparser(subparsers, cli_imagemol)
libcli.add_subparser(subparsers, cli_noSol)
libcli.add_subparser(subparsers, cli_report)
libcli.add_subparser(subparsers, cli_rmsd)
libcli.add_subparser(subparsers, cli_rmsf)
libcli.add_subparser(subparsers, cli_rot)
libcli.add_subparser(subparsers, cli_trajedit)


#ap_angle = subparsers.add_parser(
#    'angle',
#    description=cli_angle.__doc__,
#    help='Calculates the angle between a plane along the trajectory.',
#    parents=[cli_angle.ap],
#    add_help=False,
#    )
#ap_angle.set_defaults(func=cli_angle.main)

#ap_dist = subparsers.add_parser(
#    'dist',
#    description=cli_dist.__doc__,
#    help='Calculates distances between geometric centres of selections',
#    parents=[cli_dist.ap],
#    add_help=False,
#    )
#ap_dist.set_defaults(func=cli_dist.main)
#
#ap_imagemol = subparsers.add_parser(
#    'imagemol',
#    description=cli_imagemol.__doc__,
#    help='Attempts to image molecules.',
#    parents=[cli_imagemol.ap],
#    add_help=False,
#    )
#ap_imagemol.set_defaults(func=cli_imagemol.main)
#
#ap_fext = subparsers.add_parser(
#    'fext',
#    description=cli_fext.__doc__,
#    help='Extracts single frames from trajectory.',
#    parents=[cli_fext.ap],
#    add_help=False,
#    )
#ap_fext.set_defaults(func=cli_fext.main)
#
#ap_noSol = subparsers.add_parser(
#    'noSol',
#    description=cli_noSol.__doc__,
#    help='Removes solvent and extracts fist frame',
#    parents=[cli_noSol.ap],
#    add_help=False,
#    )
#ap_noSol.set_defaults(func=cli_noSol.main)
#
#ap_tedit = subparsers.add_parser(
#    'trajedit',
#    description=cli_trajedit.__doc__,
#    help='Converts traj to PDB.',
#    parents=[cli_trajedit.ap],
#    add_help=False,
#    )
#ap_tedit.set_defaults(func=cli_trajedit.main)
#
#ap_report = subparsers.add_parser(
#    'report',
#    description=cli_report.__doc__,
#    help='Reports on traj details.',
#    parents=[cli_report.ap],
#    add_help=False,
#    )
#ap_report.set_defaults(func=cli_report.main)
#
#ap_rot = subparsers.add_parser(
#    'rotations',
#    description=cli_rot.__doc__,
#    help='Calculates angular rotations across axes.',
#    parents=[cli_rot.ap],
#    add_help=False,
#    )
#ap_rot.set_defaults(func=cli_rot.main)
#
#ap_rmsd = subparsers.add_parser(
#    'rmsd',
#    description=cli_rmsd.__doc__,
#    help='Calculates and plot RMSDs',
#    parents=[cli_rmsd.ap],
#    add_help=False,
#    )
#ap_rmsd.set_defaults(func=cli_rmsd.main)
#
#ap_rmsf = subparsers.add_parser(
#    'rmsf',
#    description=cli_rmsf.__doc__,
#    help='Calculates and plot RMSFs',
#    parents=[cli_rmsf.ap],
#    add_help=False,
#    )
#ap_rmsf.set_defaults(func=cli_rmsf.main)


def _ap():
    return ap


def load_args():
    return ap.parse_args()


def maincli():
    libcli.save_command(CMDFILE, *sys.argv)

    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    args = load_args()
    log.debug(args)
    args.func(**vars(args))


if __name__ == '__main__':
    maincli()
