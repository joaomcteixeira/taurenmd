"""
q
"""
import argparse
import simtk.openmm.app as app
import mdtraj


def load_args():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        'cif_path',
        help='The .cif structure',
        )
    
    ap.add_argument(
        'traj',
        help='The trajectory',
        )
    
    ap.add_argument(
        '-o',
        '--output_pdb',
        help='Output PDB file.',
        default='production_noHOH.pdb',
        )

    ap.add_argument(
        '-d',
        '--traj_output',
        help='The trajectory output.',
        default='production_noHOH.dcd',
        )

    cmd = ap.parse_args()
    return cmd


def main():
    
    cmd = load_args()
    main_script(**vars(cmd))

def main_script(
        cif_path,
        traj,
        output_pdb='production_noHOH.pdb',
        traj_output='production_noHOH.dcd',
        ):

    mol = app.PDBxFile(cif_path)
    top = mdtraj.Topology.from_openmm(mol.topology)
    trj = mdtraj.load(traj, top=top)
    trj.remove_solvent(inplace=True)
    trj[0].save_pdb(output_pdb)
    trj.save(traj_output)

if __name__ == '__main__':
   main() 
