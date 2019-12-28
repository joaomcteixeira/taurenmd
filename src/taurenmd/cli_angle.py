"""
Plane angle calculator.

Calculates the angle of a plane against itself in the reference frame along the whole trajectory. The plane is defined by the three centres of geometry of tree given selections.
"""
import argparse

from bioplottemplates.plots import param

from taurenmd import CMDFILE, log
from taurenmd.libs import libcalc, libcli, libmda, libio
from taurenmd.logger import S, T

_help='Calculates the angle between a plane along the trajectory.'
_name='angle'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_plane_selection_arg(ap)
libcli.add_reference_frame_arg(ap)
libcli.add_slice_arg(ap)
libcli.add_plot_arg(ap)


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    libcli.save_command(CMDFILE, *sys.argv)
    main(**vars(cmd))
    return


def main(
        topology,
        trajectories,
        plane_selection,
        ref_frame=0,
        start=None,
        stop=None,
        step=None,
        plotvars=None,
        **kwargs
        ):
    log.info(T('calculating angles'))

    u = libmda.mda_load_universe(topology, *trajectories)

    frame_slice = libio.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
    log.info(S('for slice {}', frame_slice))

    log.info(T('calculating plane eq. for reference frame'))
    log.info(S('using frame: {}', ref_frame))
    u.trajectory[ref_frame]
    reference_point_1 = u.select_atoms(plane_selection[0]).center_of_geometry()
    reference_point_2 = u.select_atoms(plane_selection[1]).center_of_geometry()
    reference_point_3 = u.select_atoms(plane_selection[2]).center_of_geometry()

    ra, rb, rc, rd = libcalc.calc_plane_eq(
        reference_point_1,
        reference_point_2,
        reference_point_3,
        )
    log.info(S('the equation is {}x + {}y + {}z = {}', ra, rb, rc, rd))
    
    log.info(T('Calculating angles'))
    angles = []
    for fi, ts in zip(
            range(len(u.trajectory))[frame_slice],
            u.trajectory[frame_slice],
            ):

        point1 = u.select_atoms(plane_selection[0]).center_of_geometry()
        point2 = u.select_atoms(plane_selection[1]).center_of_geometry()
        point3 = u.select_atoms(plane_selection[2]).center_of_geometry()
        
        a, b, c, d = libcalc.calc_plane_eq(point1, point2, point3)

        try:
            angles.append(libcalc.calc_angle(ra, rb, rc, a, b, c))
        except ValueError:
            # happens when ValueError: math domain error
            # this is division by zero, means math.acos(d) is 0
            log.info(S('ValueError returned angle 0.0 found'))
            angles.append(0.0)
    
    log.info(S('calculated a total of {} angles.', len(angles)))

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
