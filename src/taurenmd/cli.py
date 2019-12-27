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


ap = libcli.CustomParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
subparsers = ap.add_subparsers(
    title='TAURENMD SUBROUTINES',
    )

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


def _ap():
    return ap


def load_args():
    return ap.parse_args()


def maincli():

    if len(sys.argv) < 2:
        ap.print_help()
        ap.exit()

    args = load_args()
    log.debug(args)

    libcli.save_command(CMDFILE, *sys.argv)
    args.func(**vars(args))


if __name__ == '__main__':
    maincli()
