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
import os
from datetime import datetime

import numpy as np

from taurenmd import _BANNER, Path
from taurenmd import core as tcore
from taurenmd import log
from taurenmd.libs import libcalc, libcli, libmda, libutil
from taurenmd.libs.libutil import make_list
from taurenmd.logger import S, T
from taurenmd.plots import labeldots


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
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
libcli.add_insort_arg(ap)
libcli.add_atom_selections_arg(ap)
libcli.add_inverted_array(ap)
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
        selections=('protein and name CA',),
        inverted_selections=None,
        export=False,
        plot=False,
        plotvars=None,
        **kwargs
        ):
    """Execute client main logic."""
    print(inverted_selections)
    log.info(T('starting RMSFs calculation'))

    topology = Path(topology)
    trajectories = [Path(t) for t in trajectories]

    u = libmda.load_universe(topology, *trajectories, insort=insort)

    frame_slice = libmda.get_frame_slices(u, start, stop, step)

    selections_list = make_list(selections)
    print(selections_list)
    if not selections_list:
        emsg = (
            'No selections defined. Check your `selections` argument input: '
            f'{selections}'
            )
        raise ValueError(emsg)

    log.info(T('calculating RMSFs'))

    labels = []
    rmsfs = []

    # Create a list of slice objects according to `inverted_selections`
    # If true, create [::-1] slice, and the data for that selection will be
    # inverted. This is useful for plotting DNA or other antiparallel chains.
    if inverted_selections:
        inv_sels = [
            slice(None, None, -1) if bool(int(inv)) else slice(None, None, None)
            for inv in inverted_selections
            ]
    else:
        inv_sels = [slice(None, None, None)] * len(selections_list)

    for i, sel in enumerate(selections_list):
        log.info(S('for sel: {}', sel))

        atom_group = u.select_atoms(sel)

        _ = libmda.draw_atom_label_from_atom_group(atom_group)[inv_sels[i]]
        labels.append(_)

        _ = libcalc.mda_rmsf(atom_group, frame_slice=frame_slice)[inv_sels[i]]
        rmsfs.append(_)

    if export:

        data_str = (
            "# Date: {1}{0}"
            "# Topology: {2}{0}"
            "# Trajectories {3}{0}"
            "# Atom,RMSF combinations{0}"
            "# {4}{0}"
            "{5}"
            ).format(
                os.linesep,
                datetime.now().strftime("%d/%m/%Y, %H:%M:%S"),
                Path(topology).resolve(),
                ','.join(f.resolve().str() for f in trajectories),
                ','.join(f'{sel},RMSF' for sel in selections_list),
                libutil.make_csv_lines_in_interleaved_manner(rmsfs, labels),
                )

        log.info(T('Saving data'))
        log.info(S('to: {}', export))
        Path(export).write_text(data_str)

    if plot:
        log.info(T("Plotting results:"))

        armsfs = np.array(rmsfs)
        ymax = np.max(armsfs)

        cli_defaults = {
            'ymax': ymax * 1.1 if ymax > 0 else ymax * 0.9,
            'filename': 'plot_rmsfs.pdf',
            'title': 'RMSFs',
            'xlabel': 'Atoms',
            'ylabel': r'RMSF ($\AA$)',
            'x_labels': [_l.split('.')[-1] for _l in labels[0]],
            'labels': selections_list,
            }

        cli_defaults.update(plotvars or dict())

        labeldots.plot(armsfs, **cli_defaults)

        log.info(S(f'saved plot: {cli_defaults["filename"]}'))

    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
