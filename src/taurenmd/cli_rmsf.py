"""
Calculates RMSFS of a selection along the trajectory.

Uses `MDAnalysis.analysis.rms.RMSF`_, read their documentation for details.

.. _MDAnalysis.analysis.rms.RMSF: https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html?highlight=rmsf#MDAnalysis.analysis.rms.RMSF
"""  # noqa: E501
import argparse
import functools
from datetime import datetime

from bioplottemplates.plots import label_dots

from taurenmd import Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libmda, libio
from taurenmd.logger import S, T

_help = 'Calculates and plots RMSFs'
_name = 'rmsf'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_atom_selection_arg(ap)
libcli.add_slice_arg(ap)
libcli.add_data_export_arg(ap)
libcli.add_plot_arg(ap)


def _ap():
    return ap


def main(
        topology,
        trajectories,
        start=None,
        stop=None,
        step=None,
        selection='all',
        plotvars=None,
        export=False,
        **kwargs
        ):
    """Execute client main logic."""
    log.info(T('starting'))
   
    u = libmda.mda_load_universe(topology, *trajectories)
   
    frame_slice = libio.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
  
    labels = []
    rmsfs = []

    atom_group = u.select_atoms(selection)
    labels = libmda.draw_atom_label_from_atom_group(atom_group)

    log.info(S('for selection: {}', selection))
    rmsfs = libcalc.mda_rmsf(
        atom_group,
        frame_slice=frame_slice,
        )

    if export:
        header = (
            "# Date: {}\n"
            "# Topology: {}\n"
            "# Trajectory : {}\n"
            "# Atom,RMSF\n"
            ).format(
                datetime.now().strftime("%d/%m/%Y, %H:%M:%S"),
                Path(topology).resolve(),
                [Path(f).resolve().str() for f in trajectory],
                )
        with open(export, 'w') as fh:
            fh.write(header)
            for label, value in zip(labels, rmsfs):
                fh.write('{},{:.5}\n'.format(label, value))

    plotvars = plotvars or dict()
    plotvars.setdefault('series_labels', selection)

    log.info(T('plot params:'))
    for k, v in plotvars.items():
        log.info(S('{} = {!r}', k, v))
    
    label_dots.plot(
        labels,
        rmsfs,
        **plotvars,
        )

    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
