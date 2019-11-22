"""
Does something.
"""
import argparse

import numpy as np

from taurenmd import log
from taurenmd.libs import libcli, libio, libmda, libutil  # noqa: F401
from taurenmd.logger import S, T


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
    '-l1',
    '--sel1',
    help='First selection.',
    default='all',
    type=str,
    )

ap.add_argument(
    '-l2',
    '--sel2',
    help='Second selection.',
    default='all',
    type=str,
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
        sel1='all',
        sel2='all',
        start=None,
        stop=None,
        step=None,
        **kwargs
        ):
    log.info(T('starting'))
    
    u = libmda.mda_load_universe(topology, *list(trajectory))

    frame_slice = libutil.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
    
    atom_sel1 = u.select_atoms(sel1)
    atom_sel2 = u.select_atoms(sel2)

    distances = np.ones(len(u.trajectory), dtype=np.float32)
    
    # https://www.mdanalysis.org/MDAnalysisTutorial/atomgroups.html
    # https://www.mdanalysis.org/docs/documentation_pages/core/groups.html#MDAnalysis.core.groups.AtomGroup.center_of_geometry
    for i, ts in enumerate(u.trajectory[frame_slice]):
        distances[i] = np.linalg.norm(
            np.subtract(
                atom_sel1.center_of_geometry(),
                atom_sel2.center_of_geometry(),
                )
            )

    print(distances)

    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()
