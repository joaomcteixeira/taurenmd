"""
Does something.
"""
import argparse
import copy
import math

import numpy as np
from bioplottemplates.plots import param


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
    '-l',
    '--plane-selection',
    help=(
        'Three selections which center of geometry defines the '
        'the plane. For example: \'segid A\' \'segidB \' \'segid C\'.'
        ' The angle will be calculated between this plane along the trajectory '
        ' to the reference frame.'
        ),
    nargs=3
    )

ap.add_argument(
    '-f',
    '--frame',
    help='Calc distances to a frame instead of relative along traj.',
    default=0,
    type=int,
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


def calc_plane_eq(p1, p2, p3):
    # http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, p3)
    return a, b, c, d


def calc_angle(a1, b1, c1, a2, b2, c2):
    # from: https://www.geeksforgeeks.org/angle-between-two-planes-in-3d/ 
    d = ( a1 * a2 + b1 * b2 + c1 * c2 )
    e1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1)
    e2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2)
    d = d / (e1 * e2)
    try:
        return math.degrees(math.acos(d))
    except ValueError:
        # happens when ValueError: math domain error
        # this is division by zero, means math.acos(d) is 0
        log.info(S('ValueError returned angle 0.0 found'))
        return 0.0


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
        plane_selection,
        frame=0,
        start=None,
        stop=None,
        step=None,
        plotvars=None,
        **kwargs
        ):
    log.info(T('starting'))

    u = libmda.mda_load_universe(topology, *list(trajectory))

    frame_slice = libutil.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
   
    log.info(T('calculating plane eq. for first frame'))
    u.trajectory[0]
    reference_point_1 = u.select_atoms(plane_selection[0]).center_of_geometry()
    reference_point_2 = u.select_atoms(plane_selection[1]).center_of_geometry()
    reference_point_3 = u.select_atoms(plane_selection[2]).center_of_geometry()

    ra, rb, rc, rd = calc_plane_eq(
        reference_point_1,
        reference_point_2,
        reference_point_3,
        )
    log.info(S('the equation is {}x + {}y + {}z = {}', ra, rb, rc, rd))
    
    log.info(T('Calculating angles'))
    angles = []
    for ts in u.trajectory[frame_slice]:
        
        point1 = u.select_atoms(plane_selection[0]).center_of_geometry()
        point2 = u.select_atoms(plane_selection[1]).center_of_geometry()
        point3 = u.select_atoms(plane_selection[2]).center_of_geometry()
        
        a, b, c, d = calc_plane_eq(point1, point2, point3)

        angles.append(calc_angle(ra, rb, rc, a, b, c))
   
    print(angles)

    if plotvars is None:
        plotvars = dict()
    
    if 'labels' not in plotvars:
        plotvars['labels'] = 'plane: {}'.format(' and '.join(plane_selection))

    log.info(T('plot params:'))
    for k, v in plotvars.items():
        log.info(S('{} = {!r}', k, v))

    param.plot(
        list(range(len(u.trajectory))[frame_slice]),
        angles,
        **plotvars,
        )

    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()
