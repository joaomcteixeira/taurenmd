"""
Does something.
"""
import argparse


from sympy import Plane, Point3D
from taurenmd import log
from taurenmd.libs import libcli, libio, libmda  # noqa: F401
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
    '--origin-selection',
    help='The selection to define the plane',
    required=True,
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
        **kwargs,
        ):
    log.info(T('starting'))
    
    u = libmda.mda_load_universe(topology, *list(trajectory))

    atom_group_origin = u.select_atoms(origin_selection)
    _3_coordinates = [sel.strip() for sel in origin_selection.split('or')]
    # p stands for point
    pA_atomG = u.select_atoms(_3_coordinates[0])
    pB_atomG = u.select_atoms(_3_coordinates[1])
    pC_atomG = u.select_atoms(_3_coordinates[2])

    u.trajectory[0]
    
    # defining the center of reference
    pABC_cog = atom_group_origin.center_of_geometry()
    log.info(T('Original Center of Geometry'))
    log.info(S('for frame: 0'))
    log.info(S('for selection: {}', origin_selection))
    log.info(S('pABC_cog: {}', pABC_cog))
    
    log.info(T('Transfering'))
    log.info(S('all coordinates of reference frame to the origin 0, 0, 0'))
    atom_group_origin.positions = atom_group_origin.positions - pABC_cog
    log.info(S('COG in origin: {}', atom_group_origin.center_of_geometry()))

    log.info(T('defining the reference axes'))
    pA_cog = pA_atomG.center_of_geometry()
    pB_cog = pB_atomG.center_of_geometry()
    pC_cog = pC_atomG.center_of_geometry()
    log.info(S('plane points definition:'))
    log.info(S('pA: {}', pA_cog))
    log.info(S('pB: {}', pB_cog))
    log.info(S('pC: {}', pC_cog))

    log.info(T('Defining Reference Plane'))
    reference_plane = Plane(
        Point3D(pA_cog),
        Point3D(pB_cog),
        Point3D(pC_cog),
        )
    log.info(S('done'))
    
    log.info(T('defining the normal vector to reference plane'))
    ref_plane_normal = reference_plane.normal_vector
    log.info(S('done'))

    log.info(T('defining the cross product vector'))
    ref_cross = Plane(
        Point3D(pA_cog),
        Point3D(atom_group_origin.center_of_geometry()),
        Point3D(ref_plane_normal),
        ).normal_vector



    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()
