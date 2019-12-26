"""
Edit a trajectory.

File format, length (that is, frames) and selection can be edited.

Uses MDAnalsysis.
"""
import argparse

import MDAnalysis as mda
from MDAnalysis.analysis import align as mdaalign

from taurenmd import Path, log
from taurenmd.libs import libio, libmda
from taurenmd.logger import S, T

_help = 'Edits trajectory in many different ways.'
_name = 'trajedit'

ap = argparse.ArgumentParser(
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
    '-i',
    '--insort',
    help=(
        'Sorts input trajectory paths according to their tail numbers, '
        'if paths are formatted: my_trajectory_#.dcd, '
        'where # is a number.'
        ),
    action='store_true',
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
    help=(
        'Atom selection for the output trajectory. '
        'Read: https://www.mdanalysis.org/docs/documentation_pages/selections.html'  # noqa: E501
        ),
    default='all',
    type=str,
    )

ap.add_argument(
    '-d',
    '--traj-output',
    help='Edited trajectory. Defaults to traj_output.xtc.',
    default='traj_output.xtc',
    type=Path,
    )

ap.add_argument(
    '-o',
    '--top-output',
    help=(
        'Topology output first frame.'
        'Defaults to --traj-output file name + _frame0.'
        ),
    default=None,
    )

ap.add_argument(
    '-O',
    '--save-frame0-topology',
    help=(
        'Oposite of -o. '
        'Do NOT save frame0 as topology. '
        'Defaults to False, that is, saves topology. '
        ),
    action='store_false',
    )

ap.add_argument(
    '-a',
    '--align',
    help='Align selection (l) to a atom subselection',
    action='store_true',
    )

ap.add_argument(
    '--align-selection',
    help=(
        'The reference atom group to which align the trajectory to. '
        'Must be subselection of --selection.'
        ),
    type=str,
    default='all',
    )

ap.add_argument(
    '-w',
    '--unwrap',
    help=(
        'Unwraps selection according to: '
        'https://www.mdanalysis.org/docs/documentation_pages/core/groups.html#MDAnalysis.core.groups.AtomGroup.unwrap'  # noqa: E501
        ),
    action='store_true',
    )

ap.add_argument(
    '--unwrap-reference',
    help=(
        'The unwrap method reference parameter. '
        'Has effect only if \'-w\' is given. '
        'Defaults to `None`.'
        ),
    default=None,
    )

ap.add_argument(
    '--unwrap-compound',
    help=(
        'The unwrap method compound parameter. '
        'Has effect only if \'-w\' is given. '
        'Defaults to `fragments`.'
        ),
    default='fragments',
    type=str,
    )

def _ap():
    return ap


def load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main(**vars(cmd))
    return


def main(
        topology,
        trajectory,
        insort=None,
        start=None,
        stop=None,
        step=None,
        selection='all',
        traj_output='traj_output.xtc',
        top_output=None,
        save_frame0_topology=True,
        unwrap=False,
        unwrap_reference=None,
        unwrap_compound='fragments',
        align=False,
        align_selection='all',
        **kwargs,
        ):
   
    log.info(T('editing trajectory'))
    
    if insort:
        trajectories = libio.sort_numbered_input(*list(trajectory))
    else:
        trajectories = list(trajectory)

    u = libmda.mda_load_universe(topology, trajectories)
    
    if unwrap:
        log.info(T('unwrapping'))
        log.info(S('set to: {}', unwrap))
        log.info(S('reference: {}', unwrap_reference))
        log.info(S('compound: {}', unwrap_compound))

    if align:
        u_top = mda.Universe(topology).select_atoms(selection)
        log.info(T('Alignment'))
        log.info(S('trajectory selection will be aligned to subselection:'))
        log.info(S('- {}', align_selection, indent=2))
    
    log.info(T('transformation'))
    log.info(S('slicing: {}::{}::{}', start, stop, step))
    sliceObj = slice(start, stop, step)

    log.info(S('selecting: {}', selection))
    selection = u.select_atoms(selection)
    log.info(S('with {} atoms', selection.n_atoms, indent=2))

    log.info(T('saving trajectory'))
    traj_output = Path(traj_output)
    log.info(S('destination: {}', traj_output.resolve().str()))

    with mda.Writer(traj_output.str(), selection.n_atoms) as W:
        for i, ts in zip(
                range(len(u.trajectory))[sliceObj],
                u.trajectory[sliceObj],
                ):
            
            log.info(S('working on frame: {}', i))
            
            if unwrap:
                log.debug(S('unwrapping', indent=2))
                selection.unwrap(
                    reference=unwrap_reference,
                    compound=unwrap_compound,
                    )

            if align:
                mdaalign.alignto(
                    selection,
                    u_top,
                    select=align_selection,
                    )

            W.write(selection)
    
    log.info(S('trajectory saved'))

    if save_frame0_topology:
        
        log.info(T('saving topology'))

        if top_output is None:
            top_output = libio.mk_frame_path(traj_output)
        else:
            top_output = Path(top_output)
        
        log.info(S('saving frame 0 to: {}', top_output.resolve()))
        with mda.Writer(Path(top_output).str(), selection.n_atoms) as W:
            for ts in u.trajectory[sliceObj][0:1]:
                if unwrap:
                    log.debug(S('unwrapping for topology', indent=2))
                    selection.unwrap(
                        reference=unwrap_reference,
                        compound=unwrap_compound,
                        )
                W.write(selection)
    
    log.info(S('Done'))
    return


if __name__ == '__main__':
    maincli()
