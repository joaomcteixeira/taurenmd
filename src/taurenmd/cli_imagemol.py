"""
Client Image Molecule with MDTraj
=================================

**Make molecules whole.**

Attempts to "Recenter and apply periodic boundary conditions to the
molecules in each frame of the trajectory." (MDTraj documentation)

**Algorithm:**

Uses `MDTraj.Trajectory.image_molecule <http://mdtraj.org/1.9.3/api/generated/mdtraj.Trajectory.html?highlight=image%20molecules#mdtraj.Trajectory.image_molecules>`_ and
`MDTraj.Topology.find_molecules <http://mdtraj.org/1.9.3/api/generated/mdtraj.Topology.html?highlight=find_molecules#mdtraj.Topology.find_molecules>`_.

1. Protocol 1:

    Performs ``mdtraj.top.find_molecules`` and ``mdtraj.traj.image_molecules``
    in the trajectory as a whole. ``anchor_molecules`` parameters
    get ``mdtraj.top.find_molecules[:1]``, and ``other_molecules``
    parameter receives ``mdtraj.top.find_molecules[1:]``.

2. Protocol 2:
    
    The same as protocol 1 but executes those steps for each frame
    separately. Frames are concatenated back to a whole trajectory
    at the end.

**Examples:**

1. Basic usage, ``-o`` saves the first frame in a separate topology file:

    >>> taurenmd imagemol top.pdb traj.dcd -d imaged.dcd -o

2. For trajectories with *non-standard* molecules you can use a TPR file.

    >>> taurenmd imagemol top.tpr traj.xtc -d imaged.xtc

3. Using protocol 2

    >>> taurenmd imagemol top.tpr traj.xtc -d imaged.xtc -i 2

**References:**

"""
import argparse
import functools

import taurenmd.core as tcore
from taurenmd import Path, log
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

_help = 'Attempts to image molecules.'
_name = 'imagemol'

ap = libcli.CustomParser(
        description=__doc__,
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


def protocol1(traj):
    """Attempts to image molecules acting on the whole traj."""
    log.info(T('running reimage protocol #1'))
    log.info(S('finding molecules'))

    mols = traj.top.find_molecules()
    log.info(S('done'))
    
    log.info(T('reimaging'))
    reimaged = traj.image_molecules(
        inplace=False,
        anchor_molecules=mols[:1],
        other_molecules=mols[1:],
        )
    log.info(S('done'))
    return reimaged


def protocol2(traj):
    """Attempts to image molecules frame by frame."""
    reimaged = []
    for frame in range(len(traj)):
        log.info(S('reimaging frame: {}', frame))
        
        mols = traj[frame].top.find_molecules()
    
        reimaged.append(
            traj[frame].image_molecules(
                inplace=False,
                anchor_molecules=mols[:1],
                other_molecules=mols[1:],
                )
            )

    log.info(S('concatenating traj frames'))
    # http://mdtraj.org/1.9.3/api/generated/mdtraj.join.html#mdtraj.join
    reimaged_traj = reimaged[0].join(reimaged[1:])

    return reimaged_traj


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

    log.info('Attempting image molecules')
    
    traj = libmdt.load_traj(topology, trajectory)
    
    protocols = {
        1: protocol1,
        2: protocol2,
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
