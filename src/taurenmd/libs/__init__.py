"""
Library package.

Taurenmd libraries, all prefixed as ``lib``, contain the functions used
by the client interfaces. Though taurenmd is designed to be used as a
command-line interface, we provide detailed documentation for its
internal functions.

Lib packages are organized by their scope, some have ``libs`` have
general scopes while other are tightly bound to the Molecular Dynamics
library with which they operate. For example:

#. ``libmda`` for MDAnalysis
#. ``libmdt`` for MDTraj

These modules store functions that relate solely to the scope of
MDAnalysis and MDTraj packages; normally these are I/O related operations.

On the other hand, we decided to organize other libraries based on the
scope of their functions regardless of the dependencies they use.
For example, ``libcalc`` contains functions calculate MD parameters,
and combine the usage of different libraries when needed.

Further instructions are provided within each module documentation.
"""
