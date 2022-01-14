"""
# Calculates RMSDs for a selection.

## Algorithm

Calculates the RMSD values along a trajectory slice for different
selections. If multiple selections are given creates a series data
for that selection.

RMSD is calculated using [libcalc.mda_rmsd](https://taurenmd.readthedocs.io/en/latest/reference/libcalc.html#taurenmd.libs.libcalc.mda_rmsd).

## Examples

Calculate RMSD of the whole system:

    taurenmd rmsd top.pdb traj.dcd -e rmsd.csv

Calculates RMSDs for different selections:

    taurenmd rmsd top.pdb traj.dcd -g 'segid A' 'segid B' -e

``-x`` exports the data to a CSV file. You can also plot the data with
the ``-v`` option:

    [...] -x rmsd.csv -v title=my-plot-title xlabel=frames ylabel=RMSDs ...

where ``[...]`` is the previous command example.

You can also use ``tmdrmsd`` instead of ``taurenmd rmsd``.

## References

"""  # noqa: E501
import argparse
import functools
from datetime import datetime

from taurenmd import _BANNER, Path
from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libcalc, libcli, libio, libmda
from taurenmd.libs.libutil import make_list
from taurenmd.logger import S, T
from taurenmd.plots import plotparams


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
    )

_help = 'Calculates RMSDs for a selection and trajectory slice.'
_name = 'rmsd'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_insort_arg(ap)
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
        insort=False,
        start=None,
        stop=None,
        step=None,
        ref_frame=0,
        selections='all',
        export=False,
        plot=False,
        plotvars=None,
        **kwargs
        ):
    """Execute main client logic."""
    log.info(T('starting'))

    u = libmda.load_universe(topology, *trajectories, insort=insort)

    frame_slice = libmda.get_frame_slices(u, start, stop, step)

    selections = make_list(selections)
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
                "# frame number,{}\n"
                ).format(
                    datetime.now(),
                    Path(topology).resolve(),
                    ', '.join(f.resolve().str() for f in trajectories),
                    ref_frame,
                    ','.join(selections),
                    ),
            )

    if plot:
        log.info(T("Plotting results:"))
        plotvars = plotvars or dict()

        xdata_in_time = plotvars.pop('xdata_in_time', 'ns')
        frame_list = libmda.get_frame_list_from_slice(u, frame_slice)
        xdata, xlabel = libmda.create_x_data(u, xdata_in_time, frame_list)

        # subtitle = 'Selections: {}'.format(' · '.join(selections))
        ymax = max(max(_r) for _r in rmsds)
        ymin = min(min(_r) for _r in rmsds)

        cli_defaults = {
            'ymax': ymax * 1.1 if ymax > 0 else ymax * 0.9,
            'ymin': ymin * 1.1 if ymin < 0 else ymin * 0.9,
            'filename': 'plot_rmsd.png',
            'title': 'RMSD',
            'xlabel': xlabel,
            'ylabel': r'RMSD ($\AA$)',
            'labels': selections,
            'dpi': 600,
            'legend': True,
            }

        cli_defaults.update(plotvars)

        plotparams.plot(xdata, rmsds, **cli_defaults)

        log.info(S(f'saved plot: {cli_defaults["filename"]}'))

    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
