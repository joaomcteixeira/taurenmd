"""
Converts a trajectory to PDB format.
"""
import argparse

import MDAnalysis as mda


ap = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

ap.add_argument(
    'trajectory',
    help='The trajectory file',
    )

ap.add_argument(
    'topology',
    help='Topology file.',
    type=str,
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
    '-o',
    '--output',
    help='Output trajectory.',
    default='traj_output.pdb',
    type=str,
    )


def load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main_script(**vars(cmd))
    return


def main(
        trajectory,
        topology,
        start=None,
        stop=None,
        step=None,
        selection='all',
        output='traj_output.pdb',
        **kwargs,
        ):
    
    print(start, stop, step, selection)

    u = mda.Universe(topology, trajectory)
    
    selection = u.select_atoms(selection)
    with mda.Writer(output, selection.n_atoms) as W:
        for ts in u.trajectory[slice(start, stop, step)]:
            W.write(selection)
    
    return


if __name__ == '__main__':
    maincli()
