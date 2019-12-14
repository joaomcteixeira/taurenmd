"""
Calculates RMSFS of a selection along the trajectory.

Uses `MDAnalysis.analysis.rms.RMSF`_, read their documentation for details.

.. _MDAnalysis.analysis.rms.RMSF: https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html?highlight=rmsf#MDAnalysis.analysis.rms.RMSF
"""
import argparse
from datetime import datetime

import numpy as np
from bioplottemplates.plots import label_dots

from taurenmd import Path, log  # noqa: F401
from taurenmd.libs import libcalc, libcli, libmda, libutil
from taurenmd.logger import S, T


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
    '--selections',
    help=(
        'The atom selection upon which calculate the RMSFs. '
        'You can give multiple selections to calculate multiple RMSFs sets. '
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


def get_atom_labels(atom_group):
    """
    Translate `atom_group` to list of representing strings
    for each atom.
    """

    labels = []
    for atom in atom_group:
        s = '{}.{}{}.{}'.format(
            atom.segment.segid,
            atom.residue.resnum,
            atom.resname,
            atom.name,
            )
        labels.append(s)

    return labels


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

    labels = []
    rmsfs = []
    for selection in selections:
        current_atom_group = u.select_atoms(selection)
        labels.extend(get_atom_labels(current_atom_group))

        log.info(S('for selection: {}', selection))
        rmsfs.extend(
            libcalc.mda_rmsf(
                current_atom_group,
                frame_slice=frame_slice,
                )
            )
   
    datatable = np.array([labels, rmsfs])#.T

    if export:
        header=(
            "# Date: {}\n'"
            "# Topology: {}\n"
            "# Trajectory : {}\n"
            "# Atom,RMSF\n"
            ).format(
                datetime.now(),
                Path(topology).resolve(),
                [Path(f).resolve().str() for f in trajectory],
                )
        with open(export, 'w') as fh:
            fh.write(header)
            for label, value in zip(labels, rmsfs):
                fh.write('{},{:.5}\n'.format(label, value))


    if plotvars is None:
        plotvars = dict()
    
    if 'labels' not in plotvars:
        plotvars['labels'] = selections

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
