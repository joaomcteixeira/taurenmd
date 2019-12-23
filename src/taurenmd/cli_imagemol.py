"""
Attempt image molecule with mdtraj.
"""

from taurenmd import Path, log
from taurenmd.libs import libcli, libio, libmdt
from taurenmd.logger import S, T

_help = 'Attempts to image molecules.'
_name = 'imagemol'

_TRAJOUTPUT = 'production_imaged.dcd'


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
    help='Output trajectory. Defaults to production_imaged.xtc.',
    default=_TRAJOUTPUT,
    )

ap.add_argument(
    '-o',
    '--top-output',
    help=(
        'File name to save the first frame of the imaged trajectory.'
        ' Defaults to the --traj-output path + \'_frame0.pdb.'
        ),
    default=None,
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


def main(
        topology,
        trajectory,
        traj_output=_TRAJOUTPUT,
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


if __name__ == '__main__':
    maincli()
