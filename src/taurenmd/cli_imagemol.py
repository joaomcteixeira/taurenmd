"""
Attempt image molecule with mdtraj.
"""

from taurenmd import Path, log
from taurenmd.libs import libcli, libio, libmdt
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
    help='Output trajectory. Defaults to production_imaged.xtc.',
    default='production_imaged.xtc',
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


def load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main(**vars(cmd))


def main(
        topology,
        trajectory,
        traj_output='production_imaged.xtc',
        top_output=None,
        **kwargs
        ):

    log.info('Attempting image molecules')
    
    trj = libmdt.mdtraj_load_traj(topology, trajectory)
   
    # use largest part as anchor
    log.info(T('finding molecules'))
    mols = trj.top.find_molecules()
    log.info(S('done'))

    log.info(T('imaging'))
    reimaged = trj.image_molecules(
        inplace=False,
        anchor_molecules=mols[:1],
        other_molecules=mols[1:],
        )
    log.info(S('done'))
    
    log.info(T('saving the output'))
    reimaged.save(traj_output)

    if top_output is None:
        top_output = libio.mk_frame_path(traj_output)
    
    reimaged[0].save(Path(top_output).with_suffix('.pdb').str())
    return


if __name__ == '__main__':
    maincli()
