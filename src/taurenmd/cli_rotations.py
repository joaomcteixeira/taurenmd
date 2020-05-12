"""
# Decompose Eurler angle rotations of a selection.

**EXPERIMENTAL PROTOCOL, RESULTS MAY NOT BE RELIABLE**

*Calculate the Roll, Pitch and Yaw angles along the trajectory.*

Read further on roll, pitch and yaw angles (Euler Angles) -
[wikipedia](https://en.wikipedia.org/wiki/Euler_angles).

Here we decompose these movements around the three different axis
centered at an origin using Quaternion rotation.

## Algorithm

Given a selection of three regions, ``selection A``, ``selection B``
and ``selection C``:

1. Centers the system to the three selections center of geometry for
every frame, this is called the *origin*,

2. Calculates a plane given by the center of geometries of the three
selections, plane ABC,

3. Defines the vector OA that goes from the origin to the center of
geometry of ``selection A``, this represents the Yaw axis.

4. Defines the normal vector to the plane ABC (ABCn), this represents the
Roll axis,

5. Defines the cross product betwee vectors OA and ABCn (AONn), this is the
Pitch axis,

6. Calculates this axis of reference for every frame

Calculating the angles:

Angles represent the right hand rotation around an axis of the sistem in
a i-frame compared to the reference frame. Quaterion distance is calculated
by [libcalc.generate_quaternion_rotations](https://taurenmd.readthedocs.io/en/latest/reference/libcalc.html#generate_quaternion_rotations) and
[libcal.sort_by_minimum_Qdistances](https://taurenmd.readthedocs.io/en/latest/reference/libcalc.html#taurenmd.libs.libcalc.sort_by_minimum_Qdistances).

### Roll

The roll angle is calculated by rotating the unitary vector OA
around vector ABCn until the Quaternion distance is the minimum
between the vector OAi (in frame) and vector OA in reference frame.

### Pitch

The pitch angle is calculated by rotating the unitary vector ABCn
around vector AONn until the Quaternion distance is the minimum
between the vector ABCni (in frame) and vector ABCn in reference frame.

### Yaw

The pitch angle is calculated by rotating the unitary vector AONn
around vector OA until the Quaternion distance is the minimum
between the vector AONni (in frame) and vector AONn in reference frame.

## Examples

In the case of an homotrimer, define the axis and the origin on the trimers:

    taurenmd rorations -z 'segid A' 'segid B' 'segid C' -x rotations.csv

## References

"""  # noqa: E501
import argparse
import functools
import pprint

import numpy as np

import taurenmd.core as tcore
from taurenmd import _BANNER, Path, log
from taurenmd.libs import libcalc, libcli, libio, libmda  # noqa: F401
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
    f'{tcore.ref_pyquaternion}'
    )

_help = 'Calculates angular rotations across axes.'
_name = 'rotations'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_plane_selection_arg(ap)
libcli.add_angle_unit_arg(ap)
libcli.add_slice_arg(ap)
libcli.add_data_export_arg(ap)


def _ap():
    return ap


def main(
        topology,
        trajectories,
        plane_selection,
        aunit='degrees',
        start=None,
        stop=None,
        step=None,
        export=False,
        **kwargs,
        ):
    """Execute main client logic."""
    log.info(T('starting'))

    topology = Path(topology)
    trajectories = [Path(t) for t in trajectories]

    u = libmda.load_universe(topology, *trajectories)

    log.info(T('transformation'))
    fSlice = libio.frame_slice(start=start, stop=stop, step=step)

    origin_selection = ' or '.join(plane_selection)
    pABC_atomG = u.select_atoms(origin_selection)
    ABC_selections = plane_selection
    # p stands for point
    pA_atomG = u.select_atoms(ABC_selections[0])
    pB_atomG = u.select_atoms(ABC_selections[1])
    pC_atomG = u.select_atoms(ABC_selections[2])

    u.trajectory[0]

    # defining the center of reference
    opABC_cog = pABC_atomG.center_of_geometry()
    log.info(T('Original Center of Geometry'))
    log.info(S('for frame: 0'))
    log.info(S('for selection: {}', origin_selection))
    log.info(S('pABC_cog: {}', opABC_cog))

    log.info(T('Transfering'))
    log.info(S('all coordinates of reference frame to the origin 0, 0, 0'))
    pABC_atomG.positions = pABC_atomG.positions - opABC_cog

    # defining the new center of reference
    pABC_cog = pABC_atomG.center_of_geometry().copy()
    log.info(S('COG in origin: {}', pABC_cog))

    log.info(T('Origin Center of Geometry'))
    log.info(S('for frame: 0'))
    log.info(S('for selection: {}', origin_selection))
    log.info(S('pABC_cog: {}', pABC_cog))


    log.info(T('defining the reference axes'))
    pA_cog = pA_atomG.center_of_geometry().copy()
    pB_cog = pB_atomG.center_of_geometry().copy()
    pC_cog = pC_atomG.center_of_geometry().copy()
    log.info(S('plane points definition:'))
    log.info(S('pA: {}', pA_cog))
    log.info(S('pB: {}', pB_cog))
    #log.info(S('pC: {}', pC_cog))

    log.info(T('defining the normal vector to reference plane'))
    ref_plane_normal = libcalc.calc_plane_normal(pA_cog, pB_cog, pABC_cog)
    log.info(S('plane normal: {}', ref_plane_normal))
    log.info(S('done'))

    log.info(T('defining the cross product vector'))
    ref_plane_cross = np.cross(pA_cog, ref_plane_normal)
    log.info(S('np cross: {}', ref_plane_cross))

    #roll_angles = []
    #pitch_angles = []
    #yaw_angles = []

    total_frames = len(u.trajectory[fSlice])
    sizet = (total_frames, 3)
    roll_vectors = np.empty(sizet, dtype=np.float64)
    pitch_vectors = np.empty(sizet, dtype=np.float64)
    yaw_vectors = np.empty(sizet, dtype=np.float64)

    for i, _ts in enumerate(u.trajectory[fSlice]):
        print(f'.. working for frame :{i}')

        pABC_cog_ts = pABC_atomG.center_of_geometry().copy()
        pABC_atomG.positions = pABC_atomG.positions - pABC_cog_ts
        #pABC_cog_ts = pABC_atomG.center_of_geometry().copy()

        ts_positions = pABC_atomG.positions.copy()

        # get roll
        pABC_atomG.positions = ts_positions + ref_plane_normal
        roll_vectors[i, :] = pA_atomG.center_of_geometry().copy()
        pABC_atomG.positions = ts_positions

        # get pitch
        pABC_atomG.positions = ts_positions + ref_plane_cross
        pitch_vectors[i, :] = libcalc.calc_plane_normal(
            pA_atomG.center_of_geometry(),
            pB_atomG.center_of_geometry(),
            ref_plane_cross,
            )
        pABC_atomG.positions = ts_positions

        # get yaw
        pABC_atomG.positions = ts_positions + pA_cog
        pA_cog_yaw_ts = pA_atomG.center_of_geometry().copy()
        normalv = libcalc.calc_plane_normal(
            pA_cog_yaw_ts,
            pB_atomG.center_of_geometry(),
            pA_cog,
            )
        yaw_vectors[i, :] = np.cross(pA_cog_yaw_ts, normalv)
        pABC_atomG.positions = ts_positions


    roll_torsion = np.degrees(libcalc.torsion_set(
        pA_cog,
        pABC_cog,
        ref_plane_normal,
        np.array([list(pC_cog)]),
        #roll_vectors,
        ))

    pitch_torsion = np.degrees(libcalc.torsion_set(
        ref_plane_normal,
        pABC_cog,
        ref_plane_cross,
        pitch_vectors,
        ))

    yaw_torsion = np.degrees(libcalc.torsion_set(
        ref_plane_cross,
        pABC_cog,
        pA_cog,
        yaw_vectors,
        ))


    if export:
        file_names = []
        for _fname in ['roll', 'pitch', 'yaw']:
            file_names.append(
                libio.add_prefix_to_path(
                    export,
                    f"{_fname}_angles_",
                    )
                )

        log.info(T('Saving data to files'))
        for data, fname in zip(
                [roll_torsion, pitch_torsion, yaw_torsion],
                file_names,
                ):

            log.info(S('saving {}', fname))
            libio.export_data_to_file(
                list(range(len(u.trajectory))[fSlice]),
                data,
                fname=fname,
                header=(
                    '# Topology: {}\n'
                    '# Trajectories: {}\n'
                    '# Plane Selection: {}\n'
                    '# frame,ange{}\n'
                    ).format(
                        topology,
                        ', '.join(t.resolve().str() for t in trajectories),
                        origin_selection,
                        aunit,
                        ),
                )
        log.info(S('done'))

    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
