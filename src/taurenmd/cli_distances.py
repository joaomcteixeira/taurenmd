"""
# Calculates distances between centers of geometry of two selections.

Distance is given in 3D XYZ coordinate space units.

## Algorithm

Distance between centers of geometry is calculated by::

    np.linalg.norm(np.subtract(coord1, coord2))

Where, ``coord*`` are the centers of geometry of each atom selection
``-l1`` and ``-l2``, respectively.
Read further on [np.linalg.norm](https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.norm.html)
and [np.subtract](https://docs.scipy.org/doc/numpy/reference/generated/numpy.subtract.html?highlight=subtract#numpy-subtract).

## Examples

Calculate the distances between two carbon alphas:

    taurenmd dist top.pdb traj.dcd -l1 'resnum 10 and name CA' -l2 'resnum 20 and name CA'

Calculate the distances between two chains:

    taurenmd dist top.pdb traj.dcd -l1 'segid A' -l2 'segid B'

``-x`` exports the data to a CSV file. You can also plot the data with
the ``-v`` option:

    [...] -x distances.csv -v title=my-plot-title xlabel=frames ylabel=degrees ...

where ``[...]`` is the previous command example.

``dist`` can be run directly as main command instead of subroutine:

    tmddist

## References

"""  # noqa: E501
import argparse
import functools

import numpy as np

from taurenmd import _BANNER, Path
from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libcli, libio, libmda
from taurenmd.libs.libutil import make_list
from taurenmd.logger import S, T
from taurenmd.plots import plotparams


__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}')

_help = 'Calculates distances between geometric centers of selections. '
_name = 'dist'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_insort_arg(ap)

ap.add_argument(
    '-l1',
    '--sel1',
    help='First selection (a single selection).',
    default='all',
    type=str,
    )

ap.add_argument(
    '-l2',
    '--sel2',
    help=(
        'Second selections. As many as desired. For example: --sel2 "segid A" '
        '"segid B and name CA".'
        ),
    nargs='*',
    default='all',
    type=str,
    )

libcli.add_slice_arg(ap)
libcli.add_data_export_arg(ap)
libcli.add_plot_arg(ap)


def _ap():
    return ap


def main(
        topology,
        trajectories,
        insort=False,
        sel1='all',
        sel2='all',
        start=None,
        stop=None,
        step=None,
        export=False,
        plot=False,
        plotvars=None,
        **kwargs
        ):
    """Execute main client logic."""
    log.info(T('measuring distances'))

    topology = Path(topology)
    trajectories = [Path(t) for t in trajectories]

    u = libmda.load_universe(topology, *trajectories, insort=insort)

    frame_slice = libmda.get_frame_slices(u, start, stop, step)

    log.info(T('defining atom seletions'))
    log.info(S('atom selection #1: {}', sel1))

    sels2 = make_list(sel2)
    log.info(S('atom selection #2: {}', ', '.join(sels2)))
    atom_sel1 = u.select_atoms(sel1)
    atom_sels2 = [u.select_atoms(_sel) for _sel in sels2]

    distances = np.empty(
        (len(u.trajectory[frame_slice]), len(sels2)),
        dtype=np.float64,
        )

    log.info(T('Calculating distances'))
    # https://www.mdanalysis.org/MDAnalysisTutorial/atomgroups.html
    # https://www.mdanalysis.org/docs/documentation_pages/core/groups.html#MDAnalysis.core.groups.AtomGroup.center_of_geometry

    with libcli.ProgressBar(distances.shape[0], suffix='frames') as pb:
        for i, _ts in enumerate(u.trajectory[frame_slice]):
            for j, _sel2 in enumerate(atom_sels2):

                distances[i, j] = np.linalg.norm(
                    np.subtract(
                        atom_sel1.center_of_geometry(),
                        atom_sels2[j].center_of_geometry(),
                        )
                    )

            pb.increment()

    log.info(S('calculated a total of {} distances.', distances.shape[0]))

    if export:
        libio.export_data_to_file(
            list(range(len(u.trajectory))[frame_slice]),
            *[distances[:, d] for d in range(distances.shape[1])],
            fname=export,
            header=(
                '# Distances between two selections centers of geomemtry\n'
                '# topology: {}\n'
                '# trajectories: {}\n'
                '# selection #1: {}\n'
                '# frame,{}\n'
                ).format(
                    topology,
                    ', '.join(t.resolve().str() for t in trajectories),
                    sel1,
                    ', '.join(sels2),
                    ),
            )

    if plot:
        log.info(T("Plotting results:"))
        plotvars = plotvars or dict()

        xdata_in_time = plotvars.pop('xdata_in_time', 'ns')
        frame_list = libmda.get_frame_list_from_slice(u, frame_slice)
        xdata, xlabel = libmda.create_x_data(u, xdata_in_time, frame_list)

        ymax = np.max(distances)
        ymin = np.min(distances)

        cli_defaults = {
            'dpi': 600,
            'filename': 'plot_distances.png',
            'labels': sels2,
            'legend': True,
            'title': f'Distance between a main selection and others\n{sel1}',
            'xlabel': xlabel,
            'ylabel': r'Distance ($\AA$)',
            'ymax': ymax * 1.1 if ymax > 0 else ymax * 0.9,
            'ymin': ymin * 1.1 if ymin < 0 else ymin * 0.9,
            }

        cli_defaults.update(plotvars)
        plotparams.plot(xdata, distances.T, **cli_defaults)

        log.info(S(f'saved plot: {cli_defaults["filename"]}'))

    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
