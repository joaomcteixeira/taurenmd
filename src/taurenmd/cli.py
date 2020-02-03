"""
taurenmd main client interface.

Usage Examples
--------------

To access to the complete list of *taurenmd commands* with a summary
information for each, execute:

    >>> taurenmd -h

or simply:

    >>> taurenmd

To see the current version number:

    >>> taurenmd -v

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

Logging
-------

*taurenmd* logs all its running activity as follows:

1. ``.taurenmd.cmd``, keeps an historic register of the taurenmd commands
run on a given folder together with a list of the research projects
that must be cited for that particular run; these are the libraries taurenmd
used to access and process the MD data. It also servers as a record for your
research project.

2. ``.taurenmd.log``, a user readable logging information, the very same that is
printed in the ``terminal`` during runtime. Overwrites previous runs.

3. ``.taurenmd.debug``, a full verbose log with all runtime information for the
LAST run. Overwrites previous runs.
"""  # noqa: E501
import argparse
import sys

# add bellow, by alphabetical order your newly developed cli
# import taurenmd.cli_NAME as cli_NAME
import taurenmd.cli_distances as cli_dist
import taurenmd.cli_fext as cli_fext
import taurenmd.cli_imagemol as cli_imagemol
import taurenmd.cli_nosol as cli_nosol
import taurenmd.cli_pangle as cli_pangle
import taurenmd.cli_report as cli_report
import taurenmd.cli_rmsd as cli_rmsd
import taurenmd.cli_rmsf as cli_rmsf
import taurenmd.cli_rotations as cli_rot
import taurenmd.cli_trajedit as cli_trajedit
from taurenmd import _INTERFACE_DESCRIPTION, log
from taurenmd.libs import libcli
from taurenmd.logger import CMDFILE


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
# add yourself to the credits
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

ap = libcli.CustomParser(
    description=_INTERFACE_DESCRIPTION + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)

subparsers = ap.add_subparsers(
    title='taurenmd subroutines',
    )

# add your client to this block following the example of the
# other clients.
# libcli.add_subparser(subparsers, cli_NAME)
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

# Done :-) There is nothing else to do here if you are implementing a new
# client interface.


def _ap():
    return ap


def load_args():
    """Load user input arguments."""
    return ap.parse_args()


def maincli():
    """Execute main client logic."""
    if len(sys.argv) < 2:
        ap.error('A subcommand is required.')

    args = load_args()
    log.debug(args)
    libcli.save_command(CMDFILE, *sys.argv)
    result = args.func(**vars(args))
    libcli.save_references()
    return result


if __name__ == '__main__':
    maincli()
