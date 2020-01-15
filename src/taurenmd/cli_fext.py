"""
# Extract trajectory frames to individual files.

Normally used to extract frames to PDB topology files so those can
be inspected independently.

## Note

    Frame number is 0-indexed.

## Examples:

Extract frames 11 to 49 (inclusive), remember frames index start at 0:

    taurenmd fext topology.pdb trajectory.dcd -s 10 -e 50

Extract the first frame:

    taurenmd fext topology.pdb trajectory.dcd -flist 0

Extract a selection of frames:

    taurenmd fext topology.pdb trajectory.dcd -flist 0,10,23,345

Frame file types can be specified:

    taurenmd fext topology.pdb trajectory.dcd -p 10 -x .dcd

  
Atom selection can be specified as well, the following extracts
only the 'segid A' atom region of the first frame. Selection rules
are as decribed for [MDAnalysis selection](https://www.mdanalysis.org/docs/documentation_pages/selections.html).

    taurenmd fext topology.pdb trajectory.xtc -flist 0 -l 'segid A'

Multiple trajectories can be given, they will be contatenated:

    taurenmd fext top.pdb traj1.xtc traj2.xtc traj3.xtc -p 10

Can also be used as main command:

    tmdfext topology.pdb ...


## References:

"""  # noqa: E501
import argparse
import functools

import taurenmd.core as tcore
from taurenmd import _BANNER, Path, log
from taurenmd.libs import libcli, libio, libmda
from taurenmd.logger import S


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
    )

_help = 'Extract trajectory frames to individual files.'
_name = 'fext'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_atom_selection_arg(ap)
libcli.add_frame_list_arg(ap)
libcli.add_slice_arg(ap)

ap.add_argument(
    '-f',
    '--prefix',
    help='String prefix for each file. Defaults to `frame_`.',
    default='frame_',
    )

ap.add_argument(
    '-x',
    '--ext',
    help='Extension of frame files. Defaulst to .pdb',
    default='.pdb',
    )


def _ap():
    return ap


def main(
        topology,
        trajectories,
        start=None,
        stop=None,
        step=None,
        flist=None,
        prefix='frame_',
        ext='pdb',
        selection='all',
        **kwargs):
    """Execute main client logic."""
    log.info('Starting...')
    
    u = libmda.load_universe(topology, *trajectories)
    
    frames_to_extract = libio.frame_list(
        len(u.trajectory),
        start=start,
        stop=stop,
        step=step,
        flist=flist,
        )

    log.info(S('extracting {} frames', len(frames_to_extract)))
    
    zeros = len(str(len(u.trajectory)))
    ext = ext.lstrip('.').strip()
    
    atom_group = u.select_atoms(selection)

    for frame in frames_to_extract:
        file_name = '{}{}.{}'.format(
            prefix,
            str(frame).zfill(zeros),
            ext,
            )

        atom_group.write(
            filename=Path(file_name),
            frames=[frame],
            )

        log.info(S('writen frame {}, to {}', frame, file_name))

    return


maincli = functools.partial(libcli.maincli, ap, main)

if __name__ == '__main__':
    maincli()
