"""
Client Calculate RMSF
=====================

**Calculates RMSFS of a selection along the trajectory slice.**

**Algorithm:**

Calculates the RMSF values along a trajectory slice for different
selections. If multiple selections are given creates a series data
for that selection.

RMSF is calculated using :py:func:`taurenmd.libs.libcalc.mda_rmsf`.

If multiple selections are given, separate calculations are performed
in sequence. Result files (data tables and plots) are exported
separately for each selection. Selections can't be overlayed easily
in a single plot because they do not share the same labels.

**Examples:**

1. Calculate RMSF of the whole system:

    >>> taurenmd rmsf top.pdb traj.dcd -e rmsf.csv

2. Calculates RMSFs for different selections:

    >>> taurenmd rmsf top.pdb traj.dcd -g 'segid A' 'segid B' -e

3. ``-x`` exports the data to a CSV file. You can also plot the data with
the ``-v`` option:

    >>> [...] -x rmsf.csv -v title=my-plot-title xlabel=frames ylabel=RMSFs ...

where ``[...]`` is the previous command example.

4. you can also use ``tmdrmsf`` instead of ``taurenmd rmsf``.

**References:**


"""  # noqa: E501
import argparse
import functools
from datetime import datetime

from bioplottemplates.plots import label_dots

from taurenmd import Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libmda, libio
from taurenmd.logger import S, T

__doc__ += (
    f'{libcli.ref_mda}'
    f'{libcli.ref_mda_selection}'
    f'{libcli.ref_plottemplates_labeldots}'
    )

_help = 'Calculates and plots RMSFs'
_name = 'rmsf'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

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
        plotvars=None,
        export=False,
        **kwargs
        ):
    """Execute client main logic."""
    log.info(T('starting RMSFs calculation'))
   
    u = libmda.mda_load_universe(topology, *trajectories)
   
    frame_slice = libio.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
  
    if selections is None:
        selections = ['all']
    
    if not isinstance(selections, list):
        raise TypeError('selections is not list type')
    
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
                fname=libio.add_prefix_to_path(
                    export,
                    f"{sel.replace(' ', '_')}_",
                    ),
                header = (
                    "# Date: {}\n"
                    "# Topology: {}\n"
                    "# Trajectories {}\n"
                    "# Atom,RMSF\n"
                    ).format(
                        datetime.now().strftime("%d/%m/%Y, %H:%M:%S"),
                        Path(topology).resolve(),
                        [Path(f).resolve().str() for f in trajectory],
                        ),
                )
    
        if plot:
            plotvars = plotvars or dict()
            plotvars.setdefault('series_labels', selection)
            plotvars['filename'] = libio.add_prefix_to_path(
                plotvars.get('filename', 'rmsf.pdf'),
                f"{sel.replace(' ', '_')}_",
                )

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
