"""
Template for command-line clients.

This file provides a template structure to build new command line clients.

Copy this file to a new file named for your command line, for example:

    >>> cp _cli_template.py cli_NAME.py

Where ``NAME`` is the name of your new client, do not use CAPS.

In your new client, replace this docstring by a proper description
of the operations performed by your client. Keep clients as simple as possible.
Two simple clients are better than a complex one. You can/should read
through the documentation of other clients for examples.

Code blocks are mandatory by default unless stated as optional.
Also, functions and variable names in the main namespace should be kept
not to break functionality with the program infrastructure.

For those optional code blocks, just delete or comment out those options
needed or not for your purpose.

Please refer to our library documentation to see all available functions.
Feel free to implement you own functions in their respective lib*.py modules.

Once it is finished, add your client to the cli.py file, instructions are
provided there as comments throughout the file.

Also, add an entry point to the setup.py file, instructions are provided
there as well, look for the 'entry_points' argument.

Test driven development is a must. You should write the test cases for
your client and newly developed functions. Visit the tests/test_clis.py,
instructions are provided by comments throughout the script. Add additional
tests on the test_lib* corresponding file if you have developed new functions.

Finally, go the ``docs/reference`` folder. There is a ``cli_template.xxx``
template file there, just copy it to a ``cli_name.rst`` file and modify
it accordingly: 1) give a title, 2) change ``NAME`` for the name of your
client. Add that file reference to the toctree in the
``docs/reference/clients.rst`` file.

You can use ``tox`` during your development to confirm code style,
documentation build, and tests. It is useful if the tests for these
three tox environments pass before you Pull Request:

    >>> tox -e check
    >>> tox -e docs
    >>> tox -e py37  # or py36 depending
"""
import argparse
import functools

import taurenmd.core as tcore
from taurenmd import log
# add here additional libraries if you need them
# libmda, libmdt, libcalc, etc...
from taurenmd.libs import libcli, libio  # noqa: F401
from taurenmd.logger import S, T


__author__ = 'Your Name'
__email__ = 'your_email@address.com'
__maintainer__ = 'Your Name'  # if you plan to give long term maintenance
# otherwise just write, 'taurenmd'.
__credits__ = ['Your Name']
__status__ = 'Production'

# This block adds to the docstring the references used by your client.
# for example, MDAnalysis or MDTraj. Remember to add a:
# **References:**
# line to your client docstring.
# taurenmd.core module contains all strings related to citation.
# please refer to that module. Bellow is an example for referencing
# MDAnalysis.
__doc__ += (
    f'{tcore.ref_mda}'
    # just keep adding as your need
    )

# updated this variables according to your client
_help = 'Short description sentence.'
_name = 'name'  # the actual command subroutine name
# this will be the name to be called with
# $ taurenmd name
# choose wisely, has to sound :-)

# this block initiates the argparse.ArgumentParser
# libcli.CustomParser is just an extension to the argparse.AP
# that provides enhanced help message. Additional information is available
# in the documentation or the libcli.py file itself.
ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

# libcli contains several functions that add common arguments to the
# argument parser, for example, output trajectory or trajectory slice
# related arguments, so you can use those functions to add those commonly
# used arguments, instead of writing them from scratch.
# Just refer to libcli for more information, those functions
# start with 'add_'.
libcli.add_topology_arg(ap)

# choose between one of the following, whether your client will/could
# receive several trajectories as input or only one. If your client
# is based on MDAnalysis you want the first choice.
# libcli.add_trajectories_arg(ap)
# libcli.add_trajectory_arg(ap)

# but, if you need to add client specific arguments, use the space bellow:
# here you just follow Python argparse interface:
#
# https://docs.python.org/3/library/argparse.html
#
# ap.add_argument(
# ...


# this is used just for documentation purposes.
def _ap():
    return ap


def main(
        topology,
        trajectory,
        # trajectories,
        # additional args...
        # ...
        **kwargs
        ):
    """Execute main client logic."""
    log.info(T('starting'))
    
    # write your logic here.
    # use the log.info and log.debug registries to log and output
    # information to the user.
    # you can always refer to other clients for examples.

    log.info(S('done'))
    return


maincli = functools.partial(libcli.maincli, ap, main)


if __name__ == '__main__':
    maincli()
