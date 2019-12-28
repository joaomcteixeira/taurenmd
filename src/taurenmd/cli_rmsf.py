"""
Calculates RMSFS of a selection along the trajectory.

Uses `MDAnalysis.analysis.rms.RMSF`_, read their documentation for details.

.. _MDAnalysis.analysis.rms.RMSF: https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html?highlight=rmsf#MDAnalysis.analysis.rms.RMSF
"""  # noqa: E501
import argparse
from datetime import datetime

from bioplottemplates.plots import label_dots

from taurenmd import CMDFILE, Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libmda, libio
from taurenmd.logger import S, T

_help = 'Calculates and plots RMSFs'
_name = 'rmsf'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_top_argument(ap)
libcli.add_traj_argument(ap)
libcli.add_slice_opt_arguments(ap)
libcli.add_export_arg(ap)
libcli.add_selection_argument(ap)
libcli.add_plot_params(ap)


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    """Main client call."""
    cmd = load_args()
    libcli.save_command(CMDFILE, *sys.argv)
    main(**vars(cmd))
    return


def main(
        topology,
        trajectory,
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
   
    u = libmda.mda_load_universe(topology, *list(trajectory))
   
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


if __name__ == '__main__':
    maincli()
