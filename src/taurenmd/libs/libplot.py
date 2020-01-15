"""Contains plotting routines.

Plotting routines can be hardcoded or linked from thirdparty projects.

References
----------

* https://python-bioplottemplates.readthedocs.io/en/latest/
"""
from bioplottemplates.plots import plotlabeld, plotparam

import taurenmd.core as tcore
from taurenmd.libs import libcli


@libcli.add_reference(tcore.ref_plottemplates_param)
def param(*args, **kwargs):
    """
    Reproduce bioplotemplates.plots.param.
    
    Visit: https://python-bioplottemplates.readthedocs.io/en/latest/
    """
    return plotparam(*args, **kwargs)


@libcli.add_reference(tcore.ref_plottemplates_labeldots)
def label_dots(*args, **kwargs):
    """
    Reproduce bioplotemplates.plots.label_dots.

    Visit: https://python-bioplottemplates.readthedocs.io/en/latest/
    """
    return plotlabeld(*args, **kwargs)
