"""
# Remove solvent from trajectory.

## Algorithm

Removes solvent from trajectory using [MDTraj.remove_solvent](http://mdtraj.org/1.9.3/api/generated/mdtraj.Trajectory.html?highlight=remove%20solvent#mdtraj.Trajectory.remove_solvent).

## Examples

Remove all solvent:

    taurenmd nosol top.pdb traj.dcd -d traj_nosol.dcd -o

Remove all solvent except for NA atoms:

    taurenmd nosol top.pdb traj.dcd -d traj_nosolNA.dcd -e Na -o

``tmdnosol`` can be used as main command:

    tmdnosol [...]

## References

"""  # noqa: E501
import argparse
import functools

import taurenmd.core as tcore
from taurenmd import _BANNER, log
from taurenmd.libs import libcli, libio, libmdt
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += tcore.ref_mdt

_help = 'Removes solvent from trajectory.'
_name = 'nosol'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectory_arg(ap)
libcli.add_traj_output_arg(ap)
libcli.add_top_output_arg(ap)

ap.add_argument(
    '-m',
    '--maintain',
    help=(
        'List of solvent residue names to maintain in trajectory. '
        'Feeds MDTraj.Trajectory.remove_solvent.exclude parameter.'
        ),
    default=None,
    nargs='+',
    )


def _ap():
    return ap


def main(
        topology,
        trajectory,
        selection=None,
        maintain=None,
        top_output='_frame0.pdb',
        traj_output='nosol.dcd',
        **kwargs
        ):
    """Execute main client logic."""
    log.info(T('Removing solvent'))

    trj = libmdt.load_traj(topology, trajectory)
    
    if selection:
        log.info(S('selecting {}', selection))
        atom_sel = trj.top.select(selection)
        trj.atom_slice(atom_sel, inplace=True)

    trj.remove_solvent(inplace=True, exclude=maintain)
    
    if top_output:
        fout = libio.parse_top_output(top_output, traj_output)
        trj[0].save(fout.str())
        log.info(S('first frame topology saved: {}', fout))
   
    log.info(T('saving trajectory'))
    log.info(S('destination: {}', traj_output))
    trj.save(traj_output)
    log.info(S('saved'))


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
