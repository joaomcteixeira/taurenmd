"""
Does something.
"""
import argparse
import copy

import numpy as np
from bioplottemplates.plots import param

from taurenmd import log
from taurenmd.libs import libcli, libio, libmda, libutil  # noqa: F401
from taurenmd.logger import S, T

_help = 'Calculates distances between geometric centres of selections. '
_name = 'dist'

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

ap.add_argument(
    '-f',
    '--frame',
    help='Calc distances to a frame instead of relative along traj.',
    default=None,
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
    '-v',
    '--plotvars',
    help=(
        'Plot variables. '
        'Example: -v xlabel=frames ylabel=RMSD color=red.'
        ),
    nargs='*',
    action=libcli.ParamsToDict,
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
        frame=None,
        plotvars=None,
        **kwargs
        ):
    log.info(T('measuring distances'))

    u = libmda.mda_load_universe(topology, *list(trajectory))

    frame_slice = libutil.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
    log.info(S('for slice {}', frame_slice))
   
    log.info(T('defining atom seletions'))
    log.info(S('atom selection #1: {}', sel1))
    log.info(S('atom selection #2: {}', sel2))
    atom_sel1 = u.select_atoms(sel1)
    atom_sel2 = u.select_atoms(sel2)

    distances = np.ones(len(u.trajectory[frame_slice]), dtype=np.float32)
    
    if frame is not None:
        u.trajectory[int(frame)]
        reference_cog = copy.deepcopy(atom_sel1.center_of_geometry())

    log.info(T('Calculating distances'))
    # https://www.mdanalysis.org/MDAnalysisTutorial/atomgroups.html
    # https://www.mdanalysis.org/docs/documentation_pages/core/groups.html#MDAnalysis.core.groups.AtomGroup.center_of_geometry
    for i, ts in enumerate(u.trajectory[frame_slice]):

        if frame is None:
            reference_cog = atom_sel1.center_of_geometry()

        distances[i] = np.linalg.norm(
            np.subtract(
                reference_cog,
                atom_sel2.center_of_geometry(),
                )
            )
    
    log.info(S('calculated a total of {} distances.', len(distances)))
    
    if plotvars is None:
        plotvars = dict()
    
    if 'labels' not in plotvars:
        plotvars['labels'] = '{} dist {}'.format(sel1, sel2)

    log.info(T('plot params:'))
    for k, v in plotvars.items():
        log.info(S('{} = {!r}', k, v))

    param.plot(
        list(range(len(u.trajectory))[frame_slice]),
        distances,
        **plotvars,
        )

    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()
