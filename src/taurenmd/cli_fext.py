"""
Extracts frames to PDB files.
"""
import argparse

from taurenmd import Path, log
from taurenmd.libs import libcli, libmda
from taurenmd.logger import S


ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    # formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    'topology',
    help='Topology file.',
    type=str,
    )

ap.add_argument(
    'trajectory',
    help=(
        'Trajectory files. If given, multiple trajectories will be'
        'contactenated by order.'
        ),
    nargs='+',
    )

ap.add_argument(
    '-s',
    '--start',
    help='Start frame for slicing.',
    default=None,
    type=int,
    )

ap.add_argument(
    '-e',
    '--stop',
    help='Stop frame for slicing: exclusive',
    default=None,
    type=int,
    )

ap.add_argument(
    '-p',
    '--step',
    help='Step value for slicing',
    default=None,
    type=int,
    )

ap.add_argument(
    '-l',
    '--flist',
    help='List of frames to extract.',
    default=False,
    nargs='+',
    )

ap.add_argument(
    '-f',
    '--prefix',
    help='String prefix for each file. Defaults to `frame_`.',
    default='frame_',
    )

ap.add_argument(
    '-x',
    '--ext',
    help='Extension of frame files. Defaulst to .pdb',
    default='.pdb',
    )


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main(**vars(cmd))
    return


def main(
        topology,
        trajectory,
        start=None,
        stop=None,
        step=None,
        flist=None,
        prefix='frame_',
        ext='pdb',
        **kwargs):
    log.info('Starting...')
    
    u = libmda.mda_load_universe(topology, *list(trajectory))
    
    lent = len(u.trajectory)

    if any((start, stop, step)):
        frames_to_extract = range(start, stop, step)

    elif flist:
        frames_to_extract = [int(i) for i in flist]
    
    else:
        frames_to_extract = range(lent)
    
    log.info(S('extracting {} frames', len(frames_to_extract)))
    
    zeros = len(str(lent))
    ext = ext.lstrip('.').strip()
    
    selection = u.select_atoms('all')

    for frame in frames_to_extract:
        file_name = '{}{}.{}'.format(
            prefix,
            str(frame).zfill(zeros),
            ext,
            )

        selection.write(
            filename=Path(file_name),
            frames=[frame],
            )

        log.info(S('writen frame {}, to {}', frame, file_name))

    return


if __name__ == '__main__':
    maincli()
