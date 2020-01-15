"""
# Calculate angular oscillation of a plane along the trajectory.

A plane is defined by the centers of geometry of three atom selection
groups. The angle between that plane in each frame and itself in the
reference frame is computed. Angle can be reported in degrees (default)
or radians.

## Algorithm

Plane equation is computed by [libcalc.calc_plane_eq](https://taurenmd.readthedocs.io/en/latest/reference/libcalc.html#taurenmd.libs.libcalc.calc_plane_eq).
Angle between planes is computed by [libcalc.calc_planes_angle](https://taurenmd.readthedocs.io/en/latest/reference/libcalc.html#taurenmd.libs.libcalc.calc_planes_angle).
Refer to our documentation page for more details.

## Examples

Given a protein of 3 subunits (chains or segids) calculate the
angle variation of a plane that crosses the protein longitudinally:

    taurenmd pangle top.pdb traj.xtc -z 'segid A' 'segid B' 'segid C' -x

``-x`` exports the data to a CSV file. You can also plot the data with
the ``-v`` option:

    [...] -v title=my-plot-title xlabel=frames ylabel=degrees ...

where ``[...]`` is the previous command example.

``pangle`` can be run directly as main command instead of subroutine:

    tmdpangle

## References

"""  # noqa: E501
import argparse
import functools

import taurenmd.core as tcore
from taurenmd import _BANNER, Path, log
from taurenmd.libs import libcalc, libcli, libio, libmda, libplot
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
    f'{tcore.ref_plottemplates_param}'
    )

_help = 'Calculate angular oscillation of a plane along the trajectory.'
_name = 'pangle'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
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
        aunit='degrees',
        ref_frame=0,
        start=None,
        stop=None,
        step=None,
        export=False,
        plot=False,
        plotvars=None,
        **kwargs
        ):
    """Execute main client logic."""
    log.info(T('calculating angles'))
    
    topology = Path(topology)
    trajectories = [Path(t) for t in trajectories]

    u = libmda.load_universe(topology, *trajectories)

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
    for _ts in u.trajectory[frame_slice]:

        point1 = u.select_atoms(plane_selection[0]).center_of_geometry()
        point2 = u.select_atoms(plane_selection[1]).center_of_geometry()
        point3 = u.select_atoms(plane_selection[2]).center_of_geometry()
        
        a, b, c, d = libcalc.calc_plane_eq(point1, point2, point3)

        angles.append(
            libcalc.calc_planes_angle(
                ra, rb, rc, a, b, c,
                aunit=aunit,
                )
            )
    
    log.info(S('calculated a total of {} angles.', len(angles)))
   
    if export:
        libio.export_data_to_file(
            list(range(len(u.trajectory))[frame_slice]),
            angles,
            fname=export,
            header=(
                '# Angular oscillation between a plane representatives\n'
                '# topology: {}\n'
                '# trajectories: {}\n'
                '# selections: {}\n'
                '# frame,angle({})\n'
                ).format(
                    topology,
                    ', '.join(t.resolve().str() for t in trajectories),
                    plane_selection,
                    aunit,
                    ),
            )

    if plot:
        plotvars = plotvars or dict()
        plotvars.setdefault(
            'labels',
            'plane: {}'.format(' and '.join(plane_selection)),
            )

        log.info(T('plot params:'))
        for k, v in plotvars.items():
            log.info(S('{} = {!r}', k, v))

        libplot.param(
            list(range(len(u.trajectory))[frame_slice]),
            angles,
            **plotvars,
            )

    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
