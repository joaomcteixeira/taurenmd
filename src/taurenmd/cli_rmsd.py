"""
Does something.
"""
import argparse
from datetime import datetime

import numpy as np
from bioplottemplates.plots import param

from taurenmd import Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libmda, libutil
from taurenmd.logger import S, T

_help = 'Calculates and plots RMSDs.'
_name = 'rmsd'

_REF_FRAME = 0


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
        'Export calculated RMSDs to CSV file. '
        'If given defaults to rmsd.csv, alternatively, '
        'you can give a specific file name.'
        ),
    default=False,
    const='rmsd.csv',
    nargs='?',
    )

ap.add_argument(
    '-l',
    '--selections',
    help=(
        'The atom selection upon which calculate the RMSDs. '
        'You can give multiple selections to calculate multiple RMSDs sets. '
        'Defauts to \'all\'.'
        ),
    type=str,
    default=None,
    nargs='+',
    )

ap.add_argument(
    '-r',
    '--ref-frame',
    help='The reference frame.',
    type=int,
    default=_REF_FRAME,
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
    cmd = load_args()
    main(**vars(cmd))
    return


def main(
        topology,
        trajectory,
        start=None,
        stop=None,
        step=None,
        ref_frame=_REF_FRAME,
        selections=None,
        plotvars=None,
        export=False,
        **kwargs
        ):
    
    log.info(T('starting'))
    
    u = libmda.mda_load_universe(topology, *list(trajectory))
    
    frame_slice = libutil.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
    
    if selections is None:
        selections = ['all']
    rmsds = []
    for selection in selections:
        rmsds.append(
            libcalc.mda_rmsd_combined_chains(
                u,
                frame_slice=frame_slice,
                selection=selection,
                ref_frame=ref_frame,
                )
            )
    if export:
        np.savetxt(
            export,
            np.array([range(len(u.trajectory))[frame_slice]] + rmsds).T,
            fmt=['%d'] + ['%.5e'] * len(rmsds),
            delimiter=',',
            newline='\n',
            header=(
                "Date: {}\n'"
                "Topology: {}\n"
                "Trajectory : {}\n"
                "ref frame: {}\n"
                "frame number, {}\n"
                ).format(
                    datetime.now(),
                    Path(topology).resolve(),
                    [Path(f).resolve().str() for f in trajectory],
                    ','.join(selections),
                    ref_frame,
                    ),
            footer='',
            comments='#',
            encoding=None,
            )

    if plotvars is None:
        plotvars = dict()
    
    if 'labels' not in plotvars:
        plotvars['labels'] = selections

    log.info(T('plot params:'))
    for k, v in plotvars.items():
        log.info(S('{} = {!r}', k, v))
    
    param.plot(
        list(range(len(u.trajectory))[frame_slice]),
        rmsds,
        **plotvars,
        )

    return


if __name__ == '__main__':
    maincli()
