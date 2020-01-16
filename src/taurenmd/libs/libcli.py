"""
Shared operations for client interfaces.

This module contains functions and classes that are shared amongst the
client interfaces. It contains also others used to enhance the user
experience.
"""
import argparse
import sys
from datetime import datetime
from functools import wraps

import taurenmd.core as tcore
from taurenmd import _BANNER, references
from taurenmd.logger import CMDFILE


def load_args(ap):
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli(ap, main):
    """
    Client main function.

    Operates when client is called directly outside the
    ``taurenmd`` client interface.
    
    - Reads input parameters
    - saves inpu command to log file
    - runs client ``main`` function
    - saves references to log file
    
    Returns
    -------
    The result value from client ``main`` function.
    """
    cmd = load_args(ap)
    save_command(CMDFILE, *sys.argv)
    result = main(**vars(cmd))
    save_references()
    return result


def add_reference(ref):
    """
    Add reference decorator.

    Example
    -------

        >>> @add_reference(str)
        >>> def myfunct():
        >>>     ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            references.add(ref)
            result = func(*args, **kwargs)
            return result
        return wrapper
    return decorator


def save_references():
    """Save used references to log file."""
    out = [f'{i}: {ref}' for i, ref in enumerate(sorted(references), start=2)]
    with open(CMDFILE, 'a') as fh:
        fh.write('References:\n')
        fh.write(f'1: {tcore.ref_taurenmd}')
        fh.write(''.join(out))


# https://stackoverflow.com/questions/4042452
class CustomParser(argparse.ArgumentParser):
    """Custom Parser class."""

    def error(self, message):
        """Present error message."""
        self.print_help()
        self.print_usage()
        sys.stderr.write('*** INPUT ERROR: {}\n'.format(message))
        sys.stderr.write('*** read the usage details above\n')
        sys.exit(2)


class ParamsToDict(argparse.Action):
    """
    Convert command-line parameters in an argument to a dictionary.

    Example
    -------
    
    Where ``-x`` is an optional argument of the command-line client
    interface.

        >>> par1=1 par2='my name' par3=[1,2,3]
        >>> {'par1': 1, 'par2': 'my name', 'par3': [1, 2, 3]}

    """

    def __call__(self, parser, namespace, values, option_string=None):
        """Execute."""
        bool_value = {
            'true': True,
            'false': False,
            }

        param_dict = {}
        for kv in values:
            # print(param_dict, kv)
            try:
                k, v = kv.split('=')
            except ValueError:
                param_dict[kv] = True
            else:
                if ',' in v:
                    vs = v.split(',')
                    try:
                        param_dict[k] = tuple(float(i) for i in vs)
                    except (ValueError, TypeError):
                        param_dict[k] = tuple(i for i in vs)

                else:
                    try:
                        param_dict[k] = float(v)
                    except (ValueError, TypeError):  # is string or list
                        param_dict[k] = bool_value.get(v.lower(), v)
            
        namespace.plotvars = param_dict
        setattr(namespace, self.dest, True)


def save_command(fname, *args):
    """
    Append the execution command to a log file.

    Parameters
    ----------
    fname : string or Path
        The file name of the log file where to append the command.

    *args : strings
        String parts that compose the command.
    """
    with open(fname, 'a') as fh:
        fh.write(
            '[{}] {}\n'.format(
                datetime.now().strftime("%d/%m/%Y, %H:%M:%S"),
                ' '.join(str(a) for a in args),
                )
            )


def add_subparser(parser, module):
    """
    Add a subcommand to a parser.

    Parameters
    ----------
    parser : `argparse.add_suparsers object <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.add_subparsers>`_
        The parser to add the subcommand to.

    module
        A python module containing the characteristics of a taurenmd
        client interface. Client interface modules require the following
        attributes: ``__doc__`` which feeds the `description argument <https://docs.python.org/3/library/argparse.html#description>`_
        of `add_parser <https://docs.python.org/3/library/argparse.html#other-utilities>`_,
        ``_help`` which feeds `help <https://docs.python.org/3/library/argparse.html#help>`_,
        ``ap`` which is an `ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_,
        and a ``main`` function, which executes the main logic of the interface.
    """  # noqa: E501
    new_ap = parser.add_parser(
        module._name,
        description=module.ap.description,
        help=module._help,
        parents=[module.ap],
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    new_ap.set_defaults(func=module.main)


# arguments list
# a: angle unit
# d: trajectory output
# e: slice stop
# g: atom selections
# i: sort input by trail int
# l: selection
# o: topology output
# p: slice step
# r: reference frame
# s: slice start
# t: framelist
# plot: plot
# x: export data to table
# z: plane selection


def add_version_arg(parser):
    """
    Add version ``-v`` option to parser.
    
    Displays a message informing the current version.
    Also accessible via ``--version``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the version argument.
    """  # noqa: E501
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        # the _BANNER contains information on the version number
        version=(
            f'{_BANNER}\n'
            f'A record of the previous versions can be found at:\n'
            f'https://taurenmd.readthedocs.io/en/latest/changelog.html\n'
            )
        )


def add_angle_unit_arg(parser):
    """
    Add angle unit selectiona argument to parser.

    Is defined by ``-a`` and ``--aunit``.

    Wether angles are to be calculated in degrees or radians.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the topology positionl argument.
    """  # noqa: E501
    parser.add_argument(
        '-a',
        '--aunit',
        help='Angular unit, either degrees or radians.',
        choices=['degrees', 'radians'],
        default='degrees',
        )


def add_insort_arg(parser):
    """
    Sort input by trail int.

    Applies :py:func:`taurenmd.libs.libio.sort_numbered_input`.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the insort argument.
    """  # noqa: E501
    parser.add_argument(
        '-i',
        '--insort',
        help=(
            'Sorts input trajectories paths according to their tail numbers, '
            'if paths are formatted as follows: my_trajectory_#.dcd, '
            'where # is a number.'
            ),
        action='store_true',
        )


def add_topology_arg(parser):
    """
    Add topology positional argument to parser.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the topology positionl argument.
    """  # noqa: E501
    parser.add_argument(
        'topology',
        help='Path to the topology file.',
        type=str,
        )


def add_trajectories_arg(parser):
    """
    Add trajectory positional argument to parser.
    
    Accepts multiple trajectory files.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the trajectory positionl argument.
    """  # noqa: E501
    parser.add_argument(
        'trajectories',
        help=(
            'Path to the trajectory files. '
            'If multiple files are given, '
            'trajectories will be concatenated by input order.'
            ),
        nargs='+',
        )


def add_trajectory_arg(parser):
    """
    Add trajectory positional argument to parser.
    
    Accepts a single trajectory file.
    
    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the trajectory positionl argument.
    """  # noqa: E501
    parser.add_argument(
        'trajectory',
        help='Path to the trajectory file.',
        )


def add_slice_arg(parser):
    """
    Add start, stop and step slicing arguments.
    
    Slicing arguments are according to `Python Slice object <https://docs.python.org/3/library/functions.html#slice>`_

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the trajectory positionl argument.
    """  # noqa: E501
    parser.add_argument(
        '-s',
        '--start',
        help=(
            'The starting index for the frame slicing. '
            'Frames are 0-indexed, so the first frame is -s 0. '
            'The starting index is inclusive. '
            'Defaults to None, considers from the beginning.'
            ),
        default=None,
        type=int,
        )

    parser.add_argument(
        '-e',
        '--stop',
        help=(
            'The ending index for the frame slicing. '
            'Frames are 0-indexed, so the last frame of a 500 frame '
            'trajectory is index 499, but because '
            'the ending index is exclusive, -e 500 is required. '
            'Defaults to None, considers to the end.'
            ),
        default=None,
        type=int,
        )

    parser.add_argument(
        '-p',
        '--step',
        help=(
            'The periodicity step value for the frame slicing, '
            '-p 10 means every 10 frames. '
            'Defaults to None, considers every 1 frame.'
            ),
        default=None,
        type=int,
        )


def add_atom_selection_arg(parser):
    """
    Add selection optional argument.

    Selection argument is a string that defines the atom selection,
    this is defined by ``-l`` and ``--selection``, and defaults to ``all``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the selection argument.
    """  # noqa: E501
    parser.add_argument(
        '-l',
        '--selection',
        help=(
            'Atom selection for the output trajectory. '
            'Selection rules are as defined by the MD analysis '
            'library used by the client interface. '
            'For instructions read the main command-line client description. '
            'Defaults to \'all\'.'
            ),
        default='all',
        type=str,
        )


def add_atom_selections_arg(parser):
    """
    Add selections optional argument.

    Selections argument is a string that defines a list of atom selections,
    this is defined by ``-g`` and ``--selections``, and defaults to ``all``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the selections argument.
    """  # noqa: E501
    parser.add_argument(
        '-g',
        '--selections',
        help=(
            'List of atom selections to operate with. '
            'Selection rules are as defined by the MD analysis '
            'library used by the client interface. '
            'For instructions read the main command-line client description. '
            'Defaults to None, uses a single selection considering all '
            'atoms. '
            "Example: -g 'segid A' 'segid B' 'name CA'"
            ),
        default=None,
        nargs='+',
        )


def add_frame_list_arg(parser):
    """
    Add frame list argument.

    Registers a list of frame numbers, is defined by ``-t`` and ``--flist``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the flist argument.
    """  # noqa: E501
    parser.add_argument(
        '-t',
        '--flist',
        help=(
            'List of frames (time steps) to consider.'
            'If applicable, this can used instead of the start, stop '
            'and step slicing arguments.'
            ),
        default=False,
        nargs='+',
        type=int,
        )


def add_plane_selection_arg(parser):
    """
    Add plane selection argument.

    Plane selection is a selection of three regions separated by 'or'
    operator.
    
    Is defined by ``-z`` and ``--plane-selection``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the plane-selection argument.
    """  # noqa: E501
    parser.add_argument(
        '-z',
        '--plane-selection',
        help=(
            'Three selection strings representing three atom regions. '
            'The plane is defined by the three centres of geometry '
            'of the three selections. For example: '
            '-z \'segid A\' \'segid B\'  \'segid C\'.'
            ),
        required=True,
        nargs=3,
        )


def add_reference_frame_arg(parser):
    """
    Add a reference frame argument.

    Reference frame is the frame to compute the parameter against.

    Depending on the client logic the reference frame might have different
    meanings.

    Is defined by ``-r`` and ``--ref-frame``.

    Defaults to ``0``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the refence-frame argument.
    """  # noqa: E501
    parser.add_argument(
        '-r',
        '--ref-frame',
        help=(
            'The frame in the trajectory that serves as '
            'reference to compute against.'
            'Defaults to 0.'
            ),
        default=0,
        type=int,
        )


def add_plot_arg(parser):
    """
    Add plot parameters.

    Plot kwargs that will be passed to the plotting function.

    Defined by ``--plot``.

    If given, plot results. Additional arguments can be given to
    specify the plot parameters.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the plot argument.
    """  # noqa: E501
    parser.add_argument(
        '--plot',
        help=(
            'Plot results. '
            'Additional arguments can be given to configure the plot '
            'style. '
            'Example: --plot xlabel=frames ylabel=RMSD color=red.'
            'Accepted plot arguments are defined by the function used '
            'to plot the result. The main description of this client '
            'which plotting function is used. '
            'Defaults to ``None``, no plot is produced.'
            ),
        nargs='*',
        default=False,
        action=ParamsToDict,
        )


def add_top_output_arg(parser):
    """
    Add argument to export first frame as topology PDB file.

    Defined by ``-o`` and ``--top-output``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the topology output argument.
    """  # noqa: E501
    parser.add_argument(
        '-o',
        '--top-output',
        help=(
            'Export edited trajectory first frame as topololy file. '
            'You can specify the exact file name, otherwise, defaults to '
            'input trajectory path + \'_frame0.pdb\'. '
            'Also, if name starts with \'_\', it is used as file suffix, '
            'if name ends with \'_\', it is used as prefix, instead.'
            ),
        default=False,
        const='_frame0.pdb',
        nargs='?',
        )


def add_traj_output_arg(parser):
    """
    Add argument to export trajectory after client modifications.

    Defined by ``-d`` and ``--traj-output``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the trajectory output argument.
    """  # noqa: E501
    parser.add_argument(
        '-d',
        '--traj-output',
        help=(
            'Modified trajectory output file name. '
            'File type will be defined by file name extension. '
            'Defaults to traj_out.dcd.'
            ),
        default='traj_out.dcd',
        )


def add_data_export_arg(parser):
    """
    Add export argument.

    Export argument configures data export to a text file in table format.

    Is defined by ``-x`` and ``--export``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the export argument.
    """  # noqa: E501
    parser.add_argument(
        '-x',
        '--export',
        help=(
            'Export calculated values to a CSV file. '
            'Defaults to \'results.csv\', alternatively, '
            'you can give a specific file name.'
            ),
        default=False,
        const='results.csv',
        nargs='?',
        )
