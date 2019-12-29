"""
Attempt image molecule with mdtraj.
"""
import argparse
import functools

from taurenmd import Path, log
from taurenmd.libs import libcli, libio, libmdt
from taurenmd.logger import S, T

_help = 'Attempts to image molecules.'
_name = 'imagemol'

ap = libcli.CustomParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

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
    """
    Attempts to image molecules acting on the whole traj.
    """
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
    """
    Attempts to image molecules frame by frame.

    .. note::

        Have not found a use for this protocol yet.

    """
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
        top_output=None,
        protocol=1,
        **kwargs
        ):

    log.info('Attempting image molecules')
    
    traj = libmdt.mdtraj_load_traj(topology, trajectory)
    
    protocols = {
        1: protocol1,
        2: protocol2,
        }

    reimaged = protocols[protocol](traj)

    log.info(T('saving the output'))
    reimaged.save_dcd(traj_output)
    log.info(S('saved trajectory: {}', traj_output))

    if top_output is None:
        top_output = libio.mk_frame_path(traj_output)
    else:
        top_output = Path(top_output)
    
    reimaged[0].save_pdb(top_output.with_suffix('.pdb').str())
    log.info(S('saving frame 0 to: {}', top_output.resolve()))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
