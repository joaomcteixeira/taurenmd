"""
Does something.
"""
import argparse
import copy
import functools

import numpy as np
from bioplottemplates.plots import param

from taurenmd import log
from taurenmd.libs import libcli, libio, libmda  # noqa: F401
from taurenmd.logger import S, T

_help = 'Calculates distances between geometric centres of selections. '
_name = 'dist'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_reference_frame_arg(ap)
libcli.add_slice_arg(ap)
libcli.add_plot_arg(ap)

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
        ref_frame=None,
        plotvars=None,
        **kwargs
        ):
    log.info(T('measuring distances'))

    u = libmda.mda_load_universe(topology, *trajectories)

    frame_slice = libio.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
    log.info(S('for slice {}', frame_slice))
   
    log.info(T('defining atom seletions'))
    log.info(S('atom selection #1: {}', sel1))
    log.info(S('atom selection #2: {}', sel2))
    atom_sel1 = u.select_atoms(sel1)
    atom_sel2 = u.select_atoms(sel2)

    distances = np.ones(len(u.trajectory[frame_slice]), dtype=np.float32)
    
    if ref_frame is not None:
        u.trajectory[int(ref_frame)]
        reference_cog = copy.deepcopy(atom_sel1.center_of_geometry())

    log.info(T('Calculating distances'))
    # https://www.mdanalysis.org/MDAnalysisTutorial/atomgroups.html
    # https://www.mdanalysis.org/docs/documentation_pages/core/groups.html#MDAnalysis.core.groups.AtomGroup.center_of_geometry
    for i, ts in enumerate(u.trajectory[frame_slice]):

        if ref_frame is None:
            reference_cog = atom_sel1.center_of_geometry()

        distances[i] = np.linalg.norm(
            np.subtract(
                reference_cog,
                atom_sel2.center_of_geometry(),
                )
            )
    
    log.info(S('calculated a total of {} distances.', len(distances)))
    
    if plotvars is None:
        plotvars = dict()
    
    if 'labels' not in plotvars:
        plotvars['labels'] = '{} dist {}'.format(sel1, sel2)

    log.info(T('plot params:'))
    for k, v in plotvars.items():
        log.info(S('{} = {!r}', k, v))

    param.plot(
        list(range(len(u.trajectory))[frame_slice]),
        distances,
        **plotvars,
        )

    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
