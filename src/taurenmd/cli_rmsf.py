"""
# Calculate RMSFS of a selection along the trajectory slice.

## Algorithm

Calculates the RMSF values along a trajectory slice for different
selections. If multiple selections are given creates a series data
for that selection.

RMSF is calculated using [libcalc.mda_rmsf](https://taurenmd.readthedocs.io/en/latest/reference/libcalc.html#taurenmd.libs.libcalc.mda_rmsf).

If multiple selections are given, separate calculations are performed
in sequence. Result files (data tables and plots) are exported
separately for each selection. Selections can't be overlayed easily
in a single plot because they do not share the same labels.

## Examples

Calculate RMSF of the whole system:

    taurenmd rmsf top.pdb traj.dcd -e rmsf.csv

Calculates RMSFs for different selections:

    taurenmd rmsf top.pdb traj.dcd -g 'segid A' 'segid B' -e

``-x`` exports the data to a CSV file. You can also plot the data with
the ``-v`` option:

    [...] -x rmsf.csv -v title=my-plot-title xlabel=frames ylabel=RMSFs ...

where ``[...]`` is the previous command example.

you can also use ``tmdrmsf`` instead of ``taurenmd rmsf``.

## References

"""  # noqa: E501
import argparse
import functools
from datetime import datetime

import taurenmd.core as tcore
from taurenmd import _BANNER, Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libio, libmda, libplot
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
    f'{tcore.ref_plottemplates_labeldots}'
    )

_help = 'Calculate RMSFS for a selection and trajectory slice.'
_name = 'rmsf'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_atom_selections_arg(ap)
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
        selections=None,
        export=False,
        plot=False,
        plotvars=None,
        **kwargs
        ):
    """Execute client main logic."""
    log.info(T('starting RMSFs calculation'))
    
    topology = Path(topology)
    trajectories = [Path(t) for t in trajectories]

    u = libmda.load_universe(topology, *trajectories)
   
    frame_slice = libio.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
  
    if selections is None:
        selections = ['all']
    
    if not isinstance(selections, list) or len(selections) == 0:
        raise TypeError('selections must be LIST with at least one element')
    
    log.info(T('calculating'))
    for sel in selections:
        labels = []
        rmsfs = []
        log.info(S('for selection: {}', sel))
        atom_group = u.select_atoms(sel)
        labels = libmda.draw_atom_label_from_atom_group(atom_group)
       
        rmsfs = libcalc.mda_rmsf(
            atom_group,
            frame_slice=frame_slice,
            )

        if export:
            libio.export_data_to_file(
                labels,
                rmsfs,
                fname=export,
                header=(
                    "# Date: {}\n"
                    "# Topology: {}\n"
                    "# Trajectories {}\n"
                    "# Atom,RMSF\n"
                    ).format(
                        datetime.now().strftime("%d/%m/%Y, %H:%M:%S"),
                        Path(topology).resolve(),
                        ', '.join(f.resolve().str() for f in trajectories),
                        ),
                )
    
        if plot:
            plotvars = plotvars or dict()
            plotvars.setdefault('series_labels', selections)

            log.info(T('plot params:'))
            for k, v in plotvars.items():
                log.info(S('{} = {!r}', k, v))
            
            libplot.label_dots(
                labels,
                rmsfs,
                **plotvars,
                )

    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
