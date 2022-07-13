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

from taurenmd import _BANNER, Path
from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libcalc, libcli, libio, libmda
from taurenmd.logger import S, T
from taurenmd.plots import plotparams


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
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
        **kwargs
        ):
    """Execute main client logic."""
    log.info(T('calculating angles'))

    topology = Path(topology)
    trajectories = [Path(t) for t in trajectories]

    u = libmda.load_universe(topology, *trajectories, insort=insort)

    frame_slice = libmda.get_frame_slices(u, start, stop, step)

    log.info(T('calculating plane eq. for reference frame'))
    log.info(S('using frame: {}', ref_frame))
    u.trajectory[ref_frame]

    log.info(T('calculating plane eq for reference points: '))
    for _sel in plane_selection:
        log.info(S(_sel))
    reference_point_1 = u.select_atoms(plane_selection[0]).center_of_geometry()
    reference_point_2 = u.select_atoms(plane_selection[1]).center_of_geometry()
    reference_point_3 = u.select_atoms(plane_selection[2]).center_of_geometry()

    ra, rb, rc, rd = libcalc.calc_plane_eq(
        reference_point_1,
        reference_point_2,
        reference_point_3,
        )
    log.info(S('the plane equation is {}x + {}y + {}z = {}', ra, rb, rc, rd))

    log.info(T('Calculating angles'))
    angles = []

    trajlen = len(u.trajectory[frame_slice])
    with libcli.ProgressBar(trajlen, suffix='frames') as pb:
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

            pb.increment()

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
        log.info(T("Plotting results:"))
        plotvars = plotvars or dict()

        xdata_in_time = plotvars.pop('xdata_in_time', 'ns')
        frame_list = libmda.get_frame_list_from_slice(u, frame_slice)
        xdata, xlabel = libmda.create_x_data(u, xdata_in_time, frame_list)

        ymax = max(angles)
        ymin = min(angles)

        cli_defaults = {
            'labels': 'plane: {}'.format(' and '.join(plane_selection)),
            'dpi': 600,
            'filename': 'plot_plane.png',
            'legend': True,
            'title': 'Plane oscilations',
            'xlabel': xlabel,
            'ylabel': aunit,
            'ymax': ymax * 1.1 if ymax > 0 else ymax * 0.9,
            'ymin': ymin * 1.1 if ymin < 0 else ymin * 0.9,
            }

        cli_defaults.update(plotvars)

        plotparams.plot(frame_list, angles, **cli_defaults)

    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
