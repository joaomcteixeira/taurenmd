"""
Client Calculate RMSDs
======================

**Calculates RMSDs of a region.**

**Algorithm:**

Calculates the RMSD values along a trajectory slice for different
selections. If multiple selections are given creates a series data
for that selection.

RMSD is calculated using :py:func:`taurenmd.libs.libcalc.mda_rmsd`.

**Examples:**

1. Calculate RMSD of the whole system:

    >>> taurenmd rmsd top.pdb traj.dcd -e rmsd.csv

2. Calculates RMSDs for different selections:

    >>> taurenmd rmsd top.pdb traj.dcd -g 'segid A' 'segid B' -e

3. ``-x`` exports the data to a CSV file. You can also plot the data with
the ``-v`` option:

    >>> [...] -x rmsd.csv -v title=my-plot-title xlabel=frames ylabel=RMSDs ...

where ``[...]`` is the previous command example.

4. you can also use ``tmdrmsd`` instead of ``taurenmd rmsd``.

**References:**


"""
import argparse
import functools
from datetime import datetime

import numpy as np
from bioplottemplates.plots import param

from taurenmd import Path, log
from taurenmd.libs import libcalc, libcli, libmda, libio
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{libcli.ref_mda}'
    f'{libcli.ref_mda_selection}'
    f'{libcli.ref_plottemplates_param}'
    )

_help = 'Calculates and plots RMSDs.'
_name = 'rmsd'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_atom_selections_arg(ap)
libcli.add_reference_frame_arg(ap)
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
        ref_frame=0,
        selections=None,
        export=False,
        plot=False,
        plotvars=None,
        **kwargs
        ):
    """Main client logic."""
    log.info(T('starting'))
    
    u = libmda.load_universe(topology, *trajectories)
    
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
        libio.export_data_to_file(
            list(range(len(u.trajectory))[frame_slice]),
            *rmsds,
            fname=export,
            delimiter=',',
            header=(
                "# Date: {}\n'"
                "# Topology: {}\n"
                "# Trajectories : {}\n"
                "# ref frame: {}\n"
                "# frame number, {}\n"
                ).format(
                    datetime.now(),
                    Path(topology).resolve(),
                    [Path(f).resolve().str() for f in trajectory],
                    ref_frame,
                    ','.join(selections),
                    ),
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


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
