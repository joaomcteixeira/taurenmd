"""
# Decompose Eurler angle rotations of a selection.

*Calculate the Roll, Pitch and Yaw angles along the trajectory.*

Read further on roll, pitch and yaw angles (Euler Angles) -
[wikipedia](https://en.wikipedia.org/wiki/Euler_angles).

Here we decompose these movements around the three different axis
centered at an origin.

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
a i-frame compared to the reference frame.

### Roll

The roll angle is given by the torsion angle between OA, origin, ABCn,
and OAi (in frame), displaced along ABCn.

### Pitch

The pitch angle is by the torsion angle between ABCn, origin, AONn, and
ABCni (in frame), displaced along AONn.

### Yaw

The yaw angle is given by the torsion angle between the AONn, origin,
OA, and AONni (in frame), displaced along OA.

## Examples

In the case of an homotrimer, define the axis and the origin on the trimers:

    taurenmd rorations -z 'segid A' 'segid B' 'segid C' -x rotations.csv

## References

"""  # noqa: E501
import argparse
import functools

import numpy as np

from taurenmd import _BANNER, Path
from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libcalc, libcli, libio, libmda, libplot  # noqa: F401
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
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
libcli.add_insort_arg(ap)
libcli.add_plane_selection_arg(ap)
libcli.add_angle_unit_arg(ap)
libcli.add_reference_frame_arg(ap)
libcli.add_slice_arg(ap)
libcli.add_data_export_arg(ap)
libcli.add_plot_arg(ap)


def _ap():
    return ap


def main(
        topology,
        trajectories,
        plane_selection,
        insort=False,
        aunit='degrees',
        ref_frame=0,
        start=None,
        stop=None,
        step=None,
        export=False,
        plot=False,
        plotvars=None,
        **kwargs,
        ):
    """Execute main client logic."""
    log.info(T('starting'))

    topology = Path(topology)
    trajectories = list(map(Path, trajectories))
    u = libmda.load_universe(topology, *trajectories, insort=False)

    log.info(T('slicing'))
    fSlice = libio.frame_slice(start=start, stop=stop, step=step)

    log.info(T('computing planes'))
    log.info(S('for reference frame: {}', ref_frame))

    ABC_selections = plane_selection  # variable renaming
    origin_selection = ' or '.join(ABC_selections)
    log.info(S('for selection: {!r}', origin_selection))

    # p stands for point
    # atomG stands for MDAnalysis AtomGroup
    pABC_atomG = u.select_atoms(origin_selection)
    pA_atomG = u.select_atoms(ABC_selections[0])
    pB_atomG = u.select_atoms(ABC_selections[1])
    pC_atomG = u.select_atoms(ABC_selections[2])

    u.trajectory[ref_frame]

    # defining the center of reference
    opABC_cog = pABC_atomG.center_of_geometry()
    log.info(S('Original pABC center of geometry: {}', opABC_cog))

    log.info(S(
        'transferring all coordinates of '
        'reference frame to the origin (0, 0, 0)'
        ))
    pABC_atomG.positions = pABC_atomG.positions - opABC_cog
    pABC_cog = pABC_atomG.center_of_geometry().copy()  # this should be zero
    log.info(S('New origin Center of Geometry'))
    log.info(S('pABC_cog: {}', pABC_cog))

    log.info(T('defining the reference axes'))
    # vector from the origin (0, 0, 0) to cog of each selection
    pA_cog = pA_atomG.center_of_geometry().copy()
    pB_cog = pB_atomG.center_of_geometry().copy()
    pC_cog = pC_atomG.center_of_geometry().copy()
    log.info(S('plane points definition:'))
    log.info(S('pA: {}', pA_cog))
    log.info(S('pB: {}', pB_cog))
    log.info(S('pC: {}', pC_cog))

    ref_plane_normal = libcalc.calc_plane_normal(pABC_cog, pA_cog, pB_cog)
    log.info(S('Normal vector to the reference plane: {}', ref_plane_normal))

    ref_plane_cross = np.cross(pA_cog, ref_plane_normal)

    log.info(S('Cross vector to pA and Normal vector: {}', ref_plane_cross))

    log.info(T('Calculating tilt angles'))
    total_frames = len(u.trajectory[fSlice])
    log.info(S('for {} frames', total_frames))
    sizet = (total_frames, 3)
    roll_vectors = np.empty(sizet, dtype=np.float64)
    pitch_vectors = np.empty(sizet, dtype=np.float64)
    yaw_vectors = np.empty(sizet, dtype=np.float64)

    with libcli.ProgressBar(total_frames, suffix='frames') as pb:
        for i, _ts in enumerate(u.trajectory[fSlice]):

            pABC_cog_ts = pABC_atomG.center_of_geometry().copy()
            pABC_atomG.positions = pABC_atomG.positions - pABC_cog_ts

            ts_positions = pABC_atomG.positions.copy()

            # get roll
            pABC_atomG.positions = ts_positions + ref_plane_normal
            roll_vectors[i, :] = pA_atomG.center_of_geometry().copy()

            # get pitch
            pABC_atomG.positions = ts_positions + ref_plane_cross
            pitch_vectors[i, :] = libcalc.calc_plane_normal(
                ref_plane_cross,
                pA_atomG.center_of_geometry(),
                pB_atomG.center_of_geometry(),
                )

            # get yaw
            pABC_atomG.positions = ts_positions + pA_cog
            pA_cog_yaw_ts = pA_atomG.center_of_geometry().copy()
            normalv = libcalc.calc_plane_normal(
                pA_cog,
                pA_cog_yaw_ts,
                pB_atomG.center_of_geometry(),
                )
            yaw_vectors[i, :] = np.cross(pA_cog_yaw_ts, normalv)
            pABC_atomG.positions = ts_positions

            pb.increment()

    log.info(S('done'))

    roll_torsion = libcalc.torsion_set(
        pA_cog,
        pABC_cog,
        ref_plane_normal,
        roll_vectors,
        )

    pitch_torsion = libcalc.torsion_set(
        ref_plane_normal,
        pABC_cog,
        ref_plane_cross,
        pitch_vectors,
        )

    yaw_torsion = libcalc.torsion_set(
        ref_plane_cross,
        pABC_cog,
        pA_cog,
        yaw_vectors,
        )

    if aunit == 'degrees':
        roll_torsion, pitch_torsion, yaw_torsion = \
            [np.degrees(a) for a in [roll_torsion, pitch_torsion, yaw_torsion]]

    if export:
        libio.export_data_to_file(
            list(range(len(u.trajectory))[fSlice]),
            roll_torsion, pitch_torsion, yaw_torsion,
            fname=export,
            header=(
                '# Topology: {}\n'
                '# Trajectories: {}\n'
                '# Plane Selection: {}\n'
                '# frame,roll,pitch,yaw [in {}]\n'
                ).format(
                    topology.resolve().str(),
                    ', '.join(t.resolve().str() for t in trajectories),
                    origin_selection,
                    aunit,
                    ),
            )

    if plot:
        plotvars = plotvars or dict()
        plotvars.setdefault('labels', ['roll', 'pitch', 'yaw'])

        libplot.param(
            list(range(len(u.trajectory))[fSlice]),
            [roll_torsion, pitch_torsion, yaw_torsion],
            **plotvars,
            )

    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
