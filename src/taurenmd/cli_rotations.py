"""
Calculate the Roll, Pitch and Yaw angles along the trajectory.

Given a selection of three regions defined as:

    selection A or selection B or selection C

calculates the roll, pitch and yaw angles of an axis of reference
calculated around that selection.

The axis of reference is calculated as follows:

    1) the centre of geometry of the selection defines the origin of the
        reference frame
       
       1.1) all frames in the trajectory are centered to that
            origin.

    2) one of the axis of the reference frame is defined by the unitary
        vectir of 'selection A'.

    3) the second axis is defined by the normal vector of the plane
        defined by the centre of geometry of the three selections
        separately.
    
    4) the last axis, is defined by the cross product of the two previous
        axis.

The above procedure is performed for each frame and the reference frame
of the first frame is stored as main reference.

Calculating the angles:

Roll)
    The roll angle is calculated by rotating the 'selection A' unitary
    vector around the 'reference normal vector' until the Quaternion
    distance is the minimum between the 'reference selection A' vector
    and the 'frame i selection A' vector.

Pitch)
    The same procedure as for Roll but the
"""
import argparse

import numpy as np

from taurenmd import log
from taurenmd.libs import libcalc, libcli, libio, libmda, libutil  # noqa: F401
from taurenmd.logger import S, T

_help = 'Calculates angular rotations across axes.'
_name = 'rotations'


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
    '--origin-selection',
    help='The selection to define the plane',
    required=True,
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
        origin_selection,
        start=None,
        stop=None,
        step=None,
        **kwargs,
        ):
    log.info(T('starting'))
    
    u = libmda.mda_load_universe(topology, *list(trajectory))

    log.info(T('transformation'))
    fSlice = libutil.frame_slice(start=start, stop=stop, step=step)
    
    pABC_atomG = u.select_atoms(origin_selection)
    ABC_selections = [sel.strip() for sel in origin_selection.split('or')]
    # p stands for point
    pA_atomG = u.select_atoms(ABC_selections[0])
    pB_atomG = u.select_atoms(ABC_selections[1])
    pC_atomG = u.select_atoms(ABC_selections[2])

    u.trajectory[0]
    
    # defining the center of reference
    pABC_cog = pABC_atomG.center_of_geometry()
    log.info(T('Original Center of Geometry'))
    log.info(S('for frame: 0'))
    log.info(S('for selection: {}', origin_selection))
    log.info(S('pABC_cog: {}', pABC_cog))
    
    log.info(T('Transfering'))
    log.info(S('all coordinates of reference frame to the origin 0, 0, 0'))
    pABC_atomG.positions = pABC_atomG.positions - pABC_cog
    log.info(S('COG in origin: {}', pABC_atomG.center_of_geometry()))

    log.info(T('defining the reference axes'))
    pA_cog = pA_atomG.center_of_geometry()
    pB_cog = pB_atomG.center_of_geometry()
    pC_cog = pC_atomG.center_of_geometry()
    log.info(S('plane points definition:'))
    log.info(S('pA: {}', pA_cog))
    log.info(S('pB: {}', pB_cog))
    log.info(S('pC: {}', pC_cog))

    log.info(T('defining the normal vector to reference plane'))
    ref_plane_normal = libcalc.calc_plane_normal(pA_cog, pB_cog, pC_cog)
    log.info(S('plane normal: {}', ref_plane_normal))
    log.info(S('done'))

    log.info(T('defining the cross product vector'))
    ref_plane_cross = np.cross(pA_cog, ref_plane_normal)
    log.info(S('np cross: {}', ref_plane_cross))
    
    roll_angles = []
    pitch_angles = []
    yaw_angles = []

    for i, ts in enumerate(u.trajectory[fSlice]):
        print(f'.. working for frame :{i}')

        pABC_cog_ts = pABC_atomG.center_of_geometry()
        pABC_atomG.positions = pABC_atomG.positions - pABC_cog_ts

        pA_cog_ts = pA_atomG.center_of_geometry()
        pB_cog_ts = pB_atomG.center_of_geometry()
        pC_cog_ts = pC_atomG.center_of_geometry()

        ts_plane_normal = libcalc.calc_plane_normal(
            pA_cog_ts,
            pB_cog_ts,
            pC_cog_ts,
            )

        ts_plane_cross = np.cross(pA_cog_ts, ts_plane_normal)

        # Calculating Quaternion Rotations
        roll_Qs_tuples = libcalc.generate_quaternion_rotations(
            ref_plane_normal,
            pA_cog_ts,
            )

        pitch_Qs_tuples = libcalc.generate_quaternion_rotations(
            ref_plane_cross,
            ts_plane_normal,
            )

        yaw_Qs_tuples = libcalc.generate_quaternion_rotations(
            pA_cog,
            ts_plane_cross,
            )

        roll_minimum = libcalc.calc_minimum_Qdistances(roll_Qs_tuples, pA_cog)
        pitch_minimum = \
            libcalc.calc_minimum_Qdistances(pitch_Qs_tuples, ref_plane_normal)
        yaw_minimum = \
            libcalc.calc_minimum_Qdistances(yaw_Qs_tuples, ref_plane_cross)
        
        roll_angles.append(round(roll_minimum.degrees, 3))
        pitch_angles.append(round(pitch_minimum.degrees, 3))
        yaw_angles.append(round(yaw_minimum.degrees, 3))

    print('... saving roll_angles.csv ...')
    libio.save_to_file(
        list(range(1, len(u.trajectory) + 1)),
        [roll_angles],
        fname='roll_angles.csv',
        )

    print('... saving pitch_angles.csv ...')
    libio.save_to_file(
        list(range(1, len(u.trajectory) + 1)),
        [pitch_angles],
        fname='pitch_angles.csv',
        )
    print('... saving yaw_angles.csv ...')
    libio.save_to_file(
        list(range(1, len(u.trajectory) + 1)),
        [yaw_angles],
        fname='yaw_angles.csv',
        )
    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()
