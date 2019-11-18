"""
Attempt image molecule with mdtraj.
"""

from taurenmd import Path, log
from taurenmd.libs import libcli, libmdt
from taurenmd.logger import S, T


ap = libcli.CustomParser()

ap.add_argument(
    'topology',
    help='The topology file',
    )

ap.add_argument(
    'trajectory',
    help='The trajectory file.',
    )

ap.add_argument(
    '-d',
    '--traj-output',
    help='Output trajectory.',
    default='production_imaged.dcd',
    )

ap.add_argument(
    '-o',
    '--top-output',
    help='First frame of the imaged trajectory.',
    default='production_imaged_frame0.pdb',
    )

ap.add_argument(
    '-p',
    '--protocol',
    help='The protocol with which reimage.',
    default=1,
    type=int,
    )


def load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main(**vars(cmd))


def protocol1(traj):
    """
    Attempts to image molecules acting on the whole traj.
    """
    log.info(T('finding molecules'))
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
    """
    reimaged = []
    for frame in range(len(traj)):
        log.info(S('reimaging frame: {}', frame))
        
        mols = traj.top.find_molecules()
    
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


def main(
        topology,
        trajectory,
        traj_output='production_imaged.dcd',
        top_output='production_imaged_frame0.pdb',
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
    reimaged.save(traj_output)
    log.info(S('saved trajectory: {}', traj_output))
    top_file = Path(top_output).with_suffix('.pdb').str()
    reimaged[0].save(top_file)
    log.info(S('saved topology for first frame: {}', top_file))
    return


if __name__ == '__main__':
    maincli()
