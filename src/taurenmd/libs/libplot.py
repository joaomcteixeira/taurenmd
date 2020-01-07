"""
Contains plotting routines.

Plotting routines can be hardcoded or linked from thirdparty projects.

References:

* https://python-bioplottemplates.readthedocs.io/en/latest/

"""
import bioplottemplates.plots as bioplots

import taurenmd.core as tcore
from taurenmd.libs import libcli


@libcli.add_reference(tcore.ref_plottemplates_param)
def param(*args, **kwargs):
    """Reproduces bioplotemplates.plots.param."""
    return bioplots.param.plot(*args, **kwargs)


@libcli.add_reference(tcore.ref_plottemplates_labeldots)
def label_dots(*args, **kwargs):
    """Reproduces bioplotemplates.plots.label_dots."""
    return bioplots.label_dots.plot(*args, **kwargs)
