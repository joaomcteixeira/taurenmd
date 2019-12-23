"""
Calculates RMSFS of a selection along the trajectory.

Uses `MDAnalysis.analysis.rms.RMSF`_, read their documentation for details.

.. _MDAnalysis.analysis.rms.RMSF: https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html?highlight=rmsf#MDAnalysis.analysis.rms.RMSF
"""  # noqa: E501
import argparse
from datetime import datetime

from bioplottemplates.plots import label_dots

from taurenmd import Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libmda, libutil
from taurenmd.logger import S, T

_help = 'Calculates and plots RMSFs'
_name = 'rmsf'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

ap.add_argument(
    'topology',
    help='Topology file.',
    type=str,
    )

ap.add_argument(
    'trajectory',
    help=(
        'Trajectory files. If given, multiple trajectories will be'
        'contactenated by order.'
        ),
    nargs='+',
    )

ap.add_argument(
    '-s',
    '--start',
    help='Start frame for slicing.',
    default=None,
    type=int,
    )

ap.add_argument(
    '-e',
    '--stop',
    help='Stop frame for slicing: exclusive',
    default=None,
    type=int,
    )

ap.add_argument(
    '-p',
    '--step',
    help='Step value for slicing',
    default=None,
    type=int,
    )

ap.add_argument(
    '-x',
    '--export',
    help=(
        'Export calculated RMSFs to CSV file. '
        'If given defaults to rmsf.csv, alternatively, '
        'you can give a specific file name.'
        ),
    default=False,
    const='rmsf.csv',
    nargs='?',
    )

ap.add_argument(
    '-l',
    '--selection',
    help='The atom selection upon which calculate the RMSFs.',
    type=str,
    default=None,
    )

ap.add_argument(
    '-v',
    '--plotvars',
    help=(
        'Plot variables. '
        'Example: -v xlabel=frames ylabel=RMSD color=red.'
        ),
    nargs='*',
    action=libcli.ParamsToDict,
    )


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    input(cmd)
    return cmd


def maincli():
    """Main client call."""
    cmd = load_args()
    main(**vars(cmd))
    return


def main(
        topology,
        trajectory,
        start=None,
        stop=None,
        step=None,
        selection=None,
        plotvars=None,
        export=False,
        **kwargs
        ):
    """Execute client main logic."""
    log.info(T('starting'))
   
    u = libmda.mda_load_universe(topology, *list(trajectory))
   
    frame_slice = libutil.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
   
    selection = selection or 'all'
  
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
