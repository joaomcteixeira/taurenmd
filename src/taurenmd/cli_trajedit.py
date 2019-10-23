"""
Edit a trajectory.

File format, length (that is, frames) and selection can be edited.

Uses MDAnalsysis.
"""
import argparse

import MDAnalysis as mda

from taurenmd import log, Path
from tauremd.logger import T, S


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
    help='The trajectory file',
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
    help='Atom selection. Read: https://www.mdanalysis.org/docs/documentation_pages/selections.html',  # noqa: E501
    default='all',
    type=str,
    )

ap.add_argument(
    '-d',
    '--traj-output',
    help='Output trajectory.',
    default='traj_output.pdb',
    type=Path,
    )

ap.add_argument(
    '-o',
    '--top-output',
    help="Output edited trajectory.",
    type=Path,
    )


def load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main_script(**vars(cmd))
    return


def main(
        topology,
        trajectory,
        start=None,
        stop=None,
        step=None,
        selection='all',
        output='traj_edited.dcd',
        **kwargs,
        ):
   
    log.info(T('editing trajectory'))

    u = libio.mda_load_universe(topology, trajectory)
    
    log.info(S(f'slicing: {start}::{stop}::{end}'))
   
    log.info(S(f'selecting: {selection}'))
    selection = u.select_atoms(selection)
    log.info(S(S(f'with {selection.n_atoms)')))

    with mda.Writer(output, selection.n_atoms) as W:
        for ts in u.trajectory[slice(start, stop, step)]:
            W.write(selection)
    
    return


if __name__ == '__main__':
    maincli()
