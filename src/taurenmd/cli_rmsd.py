"""
Does something.
"""
import argparse

from bioplottemplates.plots import param

from taurenmd import Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libio, libmda, libutil


_SELECTION = 'all'
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
    '-l',
    '--selection',
    help='The atom selection',
    type=str,
    default=_SELECTION,
    )

ap.add_argument(
    '-r',
    '--ref-frame',
    help='The refernece frame.',
    type=int,
    default=_REF_FRAME,
    )

ap.add_argument(
    '-v',
    '--plotvars',
    help='Plot variables.',
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
        selection=_SELECTION,
        plotvars=None,
        **kwargs
        ):
    
    log.info('Starting...')
    
    u = libmda.mda_load_universe(topology, *list(trajectory))
    
    frame_slice = libutil.frame_slice(
        len(u.trajectory),
        start=start,
        stop=stop,
        step=step,
        )

    rmsds_combined = libcalc.mda_rmsd_combined_chains(
        u,
        frame_slice=frame_slice,
        selection=selection,
        ref_frame=ref_frame,
        )
    
    if plotvars is None:
        plotvars = dict()

    param.plot(
        list(range(len(u.trajectory))[frame_slice]),
        rmsds_combined,
        **plotvars,
        )

    return


if __name__ == '__main__':
    maincli()
