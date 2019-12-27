"""
Shared operations for client interfaces.

This module contains functions and classes that are shared amongst the
client interfaces. It contains also others used to enhance the user
experience.
"""
import argparse
import sys
from datetime import datetime


# https://stackoverflow.com/questions/4042452
class CustomParser(argparse.ArgumentParser):
    """Custom Parser class."""
    
    def error(self, message):
        """Present error message."""
        sys.stderr.write(self.format_help())
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
            
        setattr(namespace, self.dest, param_dict)


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
                ' '.join(args),
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
        description=module.__doc__,
        help=module._help,
        parents=[module.ap],
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    new_ap.set_defaults(func=module.main)


def add_top_argument(parser):
    """
    Add topology positional argument to parser.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the topology positionl argument.
    """
    parser.add_argument(
        'topology',
        help='Path to the topology file.',
        type=str,
        )

def add_traj_argument(parser):
    """
    Add trajectory positional argument to parser.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the trajectory positionl argument.
    """
    parser.add_argument(
        'trajectory',
        help=(
            'Path to the trajectory files. '
            'If multiple files are given, '
            'trajectories will be concatenated by input order.'
            ),
        nargs='+',
        )


def add_slice_opt_arguments(parser):
    """
    Add start, stop and step slicing arguments.
    
    Slicing arguments are according to `Python Slice object <https://docs.python.org/3/library/functions.html#slice>`_

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the trajectory positionl argument.
    """
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


def add_selection_argument(parser):
    """
    Adds selection optional argument.

    Selection argument is a string that defines the atom selection,
    this is defined by ``-l`` and ``--selection``, and defaults to ``all``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the selection argument.
    """
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


def add_flist_argument(parser):
    """
    Adds frame list argument.

    Registers a list of frame numbers, is defined by ``-t`` and ``--flist``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to which add the flist argument.
    """
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
        )

