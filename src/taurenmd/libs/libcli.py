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
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
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
        )

    new_ap.set_defaults(func=module.main)
