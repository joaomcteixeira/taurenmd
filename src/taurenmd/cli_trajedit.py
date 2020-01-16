"""
# Edit a trajectory.

*In short*, takes a trajectory and apply a modification:

1. frame slicing
2. atom selection
3. unwrap
4. align
5. file format change

Saves the result to a new trajectory.

All operations are performed with `MDAnalysis <https://www.mdanalysis.org>`_.

## Examples

Changes trajectory file format

    taurenmd trajedit top.pdb traj.xtd -d traj.dcd

Extracts a part of the trajectory atoms, in this example ``segid A``,
the option ``-o`` saves the first frame of the new trajectory to a topology
file:

    taurenmd trajedit top.pdb traj.xtc -d tsegidA.dcd -o -l "segid A"

You can slice the trajectory by appending the following ``-s``, ``-e`` or
``-p`` options, this saves only every 100 frames:

    [...] -p 100

You can align the trajectory to a part of the system, for example,
align the whole system to one of its subunits:

    taurenmd trajedit top.pdb traj.dcd -d alignedA.dcd -a "segid A and name CA"

further restrain the output to a specific subselection with ``-l``:
    
    [...] -l "segid A or segid B"

``trajedit`` also implements the ``unwrap`` method from which is an
alternative approach to the ``imagemol`` client, that implements from
``MDTraj``. See references section.

    taurenmd trajedit top.pdb traj.dcd -d unwrapped.dcd -w -o unwrapped_frame0.pdb


## References

"""  # noqa: E501
import argparse
import functools

import MDAnalysis as mda
from MDAnalysis.analysis import align as mdaalign

import taurenmd.core as tcore
from taurenmd import _BANNER, Path, log
from taurenmd.libs import libcli, libio, libmda
from taurenmd.logger import S, T


__author__ = 'Joao M.C. Teixeira'
__email__ = 'joaomcteixeira@gmail.com'
__maintainer__ = 'Joao M.C. Teixeira'
__credits__ = ['Joao M.C. Teixeira']
__status__ = 'Production'

__doc__ += (
    f'{tcore.ref_mda}'
    f'{tcore.ref_mda_selection}'
    f'{tcore.ref_mda_unwrap}'
    f'{tcore.ref_mda_alignto}'
    )

_help = 'Edits trajectory in many different ways.'
_name = 'trajedit'

ap = libcli.CustomParser(
    description=_BANNER + __doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

libcli.add_version_arg(ap)
libcli.add_topology_arg(ap)
libcli.add_trajectories_arg(ap)
libcli.add_insort_arg(ap)
libcli.add_atom_selection_arg(ap)
libcli.add_slice_arg(ap)
libcli.add_traj_output_arg(ap)
libcli.add_top_output_arg(ap)


ap.add_argument(
    '-a',
    '--align',
    help='Align system to a atom group.',
    default=False,
    const='all',
    nargs='?',
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


def main(
        topology,
        trajectories,
        insort=None,
        start=None,
        stop=None,
        step=None,
        selection='all',
        traj_output='traj_out.dcd',
        top_output=None,
        unwrap=False,
        unwrap_reference=None,
        unwrap_compound='fragments',
        align=False,
        **kwargs,
        ):
    """Execute main client logic."""
    log.info(T('editing trajectory'))
    
    topology = Path(topology)
    trajectories = [Path(t) for t in trajectories]

    if insort:
        trajectories = libio.sort_numbered_input(*trajectories)

    u = libmda.load_universe(topology, *trajectories)
    
    if unwrap:
        log.info(T('unwrapping'))
        log.info(S('set to: {}', unwrap))
        log.info(S('reference: {}', unwrap_reference))
        log.info(S('compound: {}', unwrap_compound))

    if align:
        log.info(T('Alignment'))
        log.info(S('trajectory selection will be aligned to subselection:'))
        log.info(S('- {}', align, indent=2))
    
    log.info(T('transformation'))
    sliceObj = libio.frame_slice(start, stop, step)

    log.info(S('selecting: {}', selection))
    atom_selection = u.select_atoms(selection)
    log.info(S('with {} atoms', atom_selection.n_atoms, indent=2))

    log.info(T('saving trajectory'))
    traj_output = Path(traj_output)
    log.info(S('destination: {}', traj_output.resolve().str()))

    with mda.Writer(traj_output.str(), atom_selection.n_atoms) as W:
        for i, _ts in zip(
                range(len(u.trajectory))[sliceObj],
                u.trajectory[sliceObj],
                ):
            
            log.info(S('working on frame: {}', i))
            
            if unwrap:
                log.debug(S('unwrapping', indent=2))
                atom_selection.unwrap(
                    reference=unwrap_reference,
                    compound=unwrap_compound,
                    )

            if align:
                mdaalign.alignto(u, u, select=align,)

            W.write(atom_selection)
    
    log.info(S('trajectory saved'))

    if top_output:
        log.info(T('saving topology'))
        fout = libio.parse_top_output(top_output, traj_output)
        log.info(S('saving frame 0 to: {}', fout.resolve()))
        with mda.Writer(fout.str(), atom_selection.n_atoms) as W:
            for _ts in u.trajectory[sliceObj][0:1]:
                if unwrap:
                    log.debug(S('unwrapping for topology', indent=2))
                    atom_selection.unwrap(
                        reference=unwrap_reference,
                        compound=unwrap_compound,
                        )
                W.write(atom_selection)
    
    log.info(S('Done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)

if __name__ == '__main__':
    maincli()
