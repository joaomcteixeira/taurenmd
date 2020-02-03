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

import taurenmd.core as tcore
from taurenmd import _BANNER, Path, log
from taurenmd.libs import libcli, libio, libmda, libplot  # noqa: F401
from taurenmd.logger import S, T


__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
    f'{tcore.ref_plottemplates_param}'
    )

_help = 'Calculates distances between geometric centers of selections. '
_name = 'dist'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)

ap.add_argument(
    '-l1',
    '--sel1',
    help='First selection.',
    default='all',
    type=str,
    )

ap.add_argument(
    '-l2',
    '--sel2',
    help='Second selection.',
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

    u = libmda.load_universe(topology, *trajectories)
    
    frame_slice = libio.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
   
    log.info(T('defining atom seletions'))
    log.info(S('atom selection #1: {}', sel1))
    log.info(S('atom selection #2: {}', sel2))
    atom_sel1 = u.select_atoms(sel1)
    atom_sel2 = u.select_atoms(sel2)

    distances = np.ones(len(u.trajectory[frame_slice]), dtype=np.float32)
    
    log.info(T('Calculating distances'))
    # https://www.mdanalysis.org/MDAnalysisTutorial/atomgroups.html
    # https://www.mdanalysis.org/docs/documentation_pages/core/groups.html#MDAnalysis.core.groups.AtomGroup.center_of_geometry
    for i, _ts in enumerate(u.trajectory[frame_slice]):

        distances[i] = np.linalg.norm(
            np.subtract(
                atom_sel1.center_of_geometry(),
                atom_sel2.center_of_geometry(),
                )
            )
    
    log.info(S('calculated a total of {} distances.', len(distances)))
    
    if export:
        libio.export_data_to_file(
            list(range(len(u.trajectory))[frame_slice]),
            distances,
            fname=export,
            header=(
                '# Distances between two selections centers of geomemtry\n'
                '# topology: {}\n'
                '# trajectories: {}\n'
                '# selection #1: {}\n'
                '# selection #2: {}\n'
                '# frame,distance\n'
                ).format(
                    topology,
                    ', '.join(t.resolve().str() for t in trajectories),
                    sel1,
                    sel2,
                    ),
            )

    if plot:
        plotvars = plotvars or dict()
        plotvars.setdefault('labels', '{} dist {}'.format(sel1, sel2))

        log.info(T('plot params:'))
        for k, v in plotvars.items():
            log.info(S('{} = {!r}', k, v))

        libplot.param(
            list(range(len(u.trajectory))[frame_slice]),
            distances,
            **plotvars,
            )

    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
