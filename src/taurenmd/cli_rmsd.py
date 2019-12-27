"""
Does something.
"""
import argparse
from datetime import datetime

import numpy as np
from bioplottemplates.plots import param

from taurenmd import CMDFILE, Path, log
from taurenmd.libs import libcalc, libcli, libmda, libio
from taurenmd.logger import S, T

_help = 'Calculates and plots RMSDs.'
_name = 'rmsd'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_top_argument(ap)
libcli.add_traj_argument(ap)
libcli.add_slice_opt_arguments(ap)
libcli.add_selections_argument(ap)
libcli.add_export_arg(ap)
libcli.add_reference_frame(ap)
libcli.add_plot_params(ap)


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    input(cmd)
    return cmd


def maincli():
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
        ref_frame=0,
        selections=None,
        plot=False,
        plotvars=None,
        export=False,
        **kwargs
        ):
    """Main client logic."""
    log.info(T('starting'))
    
    u = libmda.load_universe(topology, *list(trajectory))
    
    frame_slice = libio.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
    
    if selections is None:
        selections = ['all']
    rmsds = []
    for selection in selections:
        rmsds.append(
            libcalc.mda_rmsd(
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
    
    if plot:

        plotvars = plotvars or dict()
        plotvars.setdefault('labels', selections)

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
