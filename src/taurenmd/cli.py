"""
Usage Examples
==============
To access to the complete list of *taurenmd commands* with a summary
information for each, execute:

    >>> taurenmd

Using ``trajedit`` as an example, lets inspect its functionality:

    >>> taurenmd trajedit -h

With ``trajedit`` you can edit a trajectory in many different ways.
For example, convert a trajectory to another format:

    >>> taurenmd trajedit topology.pdb trajectory.xtc -d new_trajectory.dcd

The above command reads the original ``trajectory.xtc`` file and outputs
the new ``new_trajectory.dcd``. You can also use ``trajedit`` to reduce
the trajectory size, say by slicing every 10 frames:

    >>> taurenmd trajedit topology.pdb trajetory.xtc -d traj_p10.xtc -p 10

the ``-p`` option refers to the slicing step size, in this case ``10`` -
reads every 10 frames. Likewise, you can pass a *start* (``-s``) and an
*end* (``-e``) arguments:

    >>> taurenmd trajedit topology.pdb trajectory.xtc -d traj_s50_e500_p10.xtc -s 50 -e 500 -p 10

Also, you can extract an Atom Selection from a trajectory to a new trajectory
file. The example bellow creates a new trajectory from the input one containing
only atoms belonging to chain A. In cases like this it is useful to extract
the atom selection as an independent topology file.

    >>> taurenmd trajedit top.pdb traj.xtc -l 'segid A' -d chainA.xtc -o chainA_topology.pdb

You can also use ``trajedit`` to extract a specific frame from a trajectory:

    >>> taurenmd trajedit topology.pdb trajectory.xtc -d frame40.pdb -s 40 -e 41

but, for this example, you could instead use the ``fext`` interface:

    >>> taurenmd fext topology.pdb trajectory.xtc -f 40 -x .pdb -f frame_

Each an every ``taurenmd`` sub command is available directly as a main
routine by prefixing a ``tmd`` to its name, for example:

    >>> taurenmd trajedit
    >>> # equals to
    >>> tmdtrajedit
"""
# link to logo
# http://patorjk.com/software/taag/#p=display&h=0&f=Epic&t=taurenmd
import argparse
import sys

import taurenmd.cli_pangle as cli_pangle
import taurenmd.cli_distances as cli_dist
import taurenmd.cli_fext as cli_fext
import taurenmd.cli_imagemol as cli_imagemol
import taurenmd.cli_nosol as cli_nosol
import taurenmd.cli_report as cli_report
import taurenmd.cli_rmsd as cli_rmsd
import taurenmd.cli_rmsf as cli_rmsf
import taurenmd.cli_rotations as cli_rot
import taurenmd.cli_trajedit as cli_trajedit
from taurenmd import __version__, _BANNER, _DOCSTRING, log
from taurenmd.libs import libcli
from taurenmd.logger import CMDFILE


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-v',
    '--version',
    action='version',
    version=(
        f'{_BANNER}\n'
        f'version: {__version__}\n\n'
        f'to see the list of all versions visit: https://github.com/joaomcteixeira/taurenmd/releases\n'
        )
    )

subparsers = ap.add_subparsers(
    title='taurenmd subroutines',
    )

libcli.add_subparser(subparsers, cli_pangle)
libcli.add_subparser(subparsers, cli_dist)
libcli.add_subparser(subparsers, cli_fext)
libcli.add_subparser(subparsers, cli_imagemol)
libcli.add_subparser(subparsers, cli_nosol)
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
        sys.stdout.write(_BANNER)
        sys.stdout.write(_DOCSTRING)
        sys.stdout.write('\n')
        ap.print_help()
        ap.print_usage()
        ap.exit()

    args = load_args()
    log.debug(args)

    libcli.save_command(CMDFILE, *sys.argv)
    args.func(**vars(args))


if __name__ == '__main__':
    maincli()
