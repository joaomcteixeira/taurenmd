"""
Copies the trajectory without solvent and extracts the first frame.
"""
import argparse

from taurenmd import CMDFILE
from taurenmd.libs import libcli, libmdt

_help = 'Removes solvent and extracts first frame'
_name = 'noSol'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_top_argument(ap)
libcli.add_single_traj_argument(ap)
libcli.add_trajout_arg(ap)
libcli.add_topout_arg(ap)


def load_args():
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    libcli.save_command(CMDFILE, *sys.argv)
    main(**vars(cmd))


def main(
        topology,
        trajectory,
        selection=None,
        output_pdb='production_noSol.pdb',
        traj_output='production_noSol.dcd',
        **kwargs
        ):

    trj = libmdt.mdtraj_load_traj(topology, trajectory)
    
    if selection:
        atom_sel = trj.top.select(selection)
        trj.atom_slice(atom_sel, inplace=True)

    trj.remove_solvent(inplace=True)
    trj[0].save_pdb(output_pdb)
    trj.save(traj_output)


if __name__ == '__main__':
    maincli()
