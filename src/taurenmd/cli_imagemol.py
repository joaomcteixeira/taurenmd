"""
# Make molecules whole.

Attempts to "Recenter and apply periodic boundary conditions to the
molecules in each frame of the trajectory." (MDTraj documentation)

## Algorithm

Uses [MDTraj.Trajectory.image_molecule](http://mdtraj.org/1.9.3/api/generated/mdtraj.Trajectory.html?highlight=image%20molecules#mdtraj.Trajectory.image_molecules) and
[MDTraj.Topology.find_molecules](http://mdtraj.org/1.9.3/api/generated/mdtraj.Topology.html?highlight=find_molecules#mdtraj.Topology.find_molecules).

### Protocol 1

Performs ``mdtraj.top.find_molecules`` and ``mdtraj.traj.image_molecules``
in the trajectory as a whole. ``anchor_molecules`` parameters
get ``mdtraj.top.find_molecules[:1]``, and ``other_molecules``
parameter receives ``mdtraj.top.find_molecules[1:]``.

### Protocol 2
    
The same as protocol 1 but executes those steps for each frame
separately. Frames are concatenated back to a whole trajectory
at the end.

## Examples

Basic usage, ``-o`` saves the first frame in a separate topology file:

    taurenmd imagemol top.pdb traj.dcd -d imaged.dcd -o

For trajectories with *non-standard* molecules you can use a TPR file.

    taurenmd imagemol top.tpr traj.xtc -d imaged.xtc

Using protocol 2

    taurenmd imagemol top.tpr traj.xtc -d imaged.xtc -i 2

## References

"""  # noqa: E501
import argparse
import functools

import taurenmd.core as tcore
from taurenmd import _BANNER, log
from taurenmd.libs import libcli, libio, libmdt
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mdt}'
    )

_help = 'Attempt to image molecules.'
_name = 'imagemol'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectory_arg(ap)
libcli.add_traj_output_arg(ap)
libcli.add_top_output_arg(ap)

ap.add_argument(
    '-i',
    '--protocol',
    help=(
        'The protocol with which reimage. '
        'Read main command description for details.'
        ),
    default=1,
    type=int,
    )


def _ap():
    return ap


def main(
        topology,
        trajectory,
        traj_output='imaged.dcd',
        top_output=False,
        protocol=1,
        **kwargs
        ):
    """Execute main client logic."""
    log.info('Attempting image molecules')
    
    traj = libmdt.load_traj(topology, trajectory)
    
    protocols = {
        1: libmdt.imagemol_protocol1,
        2: libmdt.imagemol_protocol2,
        }

    reimaged = protocols[protocol](traj)

    log.info(T('saving the output'))
    reimaged.save(traj_output)
    log.info(S('saved trajectory: {}', traj_output))

    if top_output:
        fout = libio.parse_top_output(top_output, traj_output)
        reimaged[0].save(fout.str())
        log.info(S('saving frame 0 to: {}', fout.resolve()))

    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
