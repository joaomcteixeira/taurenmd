"""
Client Plane Angular Oscillations
=================================

**Calculate angular oscillation of a plane along the trajectory.**

A plane is defined by the centers of geometry of three atom selection
groups. The angle between that plane in each frame and itself in the
reference frame is computed. Angle can be reported in degrees (default)
or radians.

**Algorithm:**

Plane equation is computed by :py:func:taurenmnd.libs.libcalc.calc_plane_eq.
Angle between planes is computed by :py:func:libcalc.calc_planes_angle.
Refer to our documentation page for more details.

**Examples:**

1. Given a protein of 3 subunits (chains or segids) calculate the
angle variation of a plane that crosses the protein longitudinally:

    >>> taurenmd pangle top.pdb traj.xtc -z 'segid A' 'segid B' 'segid C' -x

2. ``-x`` exports the data to a CSV file. You can also plot the data with
the ``-v`` option:

    >>> [...] -v title=my-plot-title xlabel=frames ylabel=degrees ...

where ``[...]`` is the previous command example.

3. ``pangle`` can be run directly as main command instead of subroutine:

    >>> tmdpangle

**References:**

* MD data is accessed using `MDAnalysis <https://www.mdanalysis.org>`_. Therefore, selection commands follow MDAnalysis `selection nomenclature <https://www.mdanalysis.org/docs/documentation_pages/selections.html#>`_.
* plotting is performed by `python-bioplottemplates plot param function <https://python-bioplottemplates.readthedocs.io/en/latest/reference/plots.html#bioplottemplates.plots.param.plot>`_. 
"""
import argparse
import functools

from bioplottemplates.plots import param

from taurenmd import log
from taurenmd.libs import libcalc, libcli, libmda, libio
from taurenmd.logger import S, T

_help='Calculates the angle between a plane along the trajectory.'
_name='pangle'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_plane_selection_arg(ap)
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
        plot=False,
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
            angles.append(
                libcalc.calc_planes_angle(
                    ra, rb, rc, a, b, c,
                    aunit=aunit,
                    )
                )
        except ValueError:
            # happens when ValueError: math domain error
            # this is division by zero, means math.acos(d) is 0
            log.info(S('ValueError returned angle 0.0 found'))
            angles.append(0.0)
    
    log.info(S('calculated a total of {} angles.', len(angles)))
   
    if export:
        libio.export_data_to_file(
            list(range(len(u.trajectory))[frame_slice]),
            angles,
            fname=export,
            header=(
                '# Angular oscillation between a plane representatives\n'
                f'# topology: {topology}\n',
                f'# trajectories: {", ".join(trajectories)}\n'
                f'# selections: {plane_selection}\n'
                f'# frame,angle({aunit})\n'
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

        param.plot(
            list(range(len(u.trajectory))[frame_slice]),
            angles,
            **plotvars,
            )

    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
