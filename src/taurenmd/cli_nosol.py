"""
Copies the trajectory without solvent and extracts the first frame.
"""
import argparse
import functools

from taurenmd.libs import libcli, libmdt

_help = 'Removes solvent and extracts first frame'
_name = 'nosol'

ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_topology_arg(ap)
libcli.add_trajectory_arg(ap)
libcli.add_traj_output_arg(ap)
libcli.add_top_output_arg(ap)


def _ap():
    return ap


def main(
        topology,
        trajectory,
        selection=None,
        output_pdb='noSol_frame0.pdb',
        traj_output='noSol.dcd',
        **kwargs
        ):

    trj = libmdt.mdtraj_load_traj(topology, trajectory)
    
    if selection:
        atom_sel = trj.top.select(selection)
        trj.atom_slice(atom_sel, inplace=True)

    trj.remove_solvent(inplace=True)
    trj[0].save_pdb(output_pdb)
    trj.save(traj_output)


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
