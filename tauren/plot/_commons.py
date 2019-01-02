"""
Common functions that serve plotting routines.

Copyright © 2018-2019 Tauren-MD Project

Contributors to this file:
- João M.C. Teixeira (https://github.com/joaomcteixeira)

Tauren-MD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Tauren-MD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tauren-MD. If not, see <http://www.gnu.org/licenses/>.
"""

import itertools as it
from matplotlib import colors as mcolors

from tauren import tlog
from tauren import tcommons

log = tlog.get_log(__name__)


def _fig_size(nsubplots, ncols=1, irow=4.8, icol=6.4):
    """
    Return the figure size (width, height) based on the number of
    subplots, number of desired columns (optional, default: 1) and
    inches per row and inches per col.
    """
    
    # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.figure.html
    
    width = ncols * icol
    height = (nsubplots // ncols) * irow
    
    fs = (width, height)
    
    log.debug(f"returning: {fs}")
    
    return fs


def _get_all_mcolors():
    """
    Returns a dictionary of colors where keys are names and values
    HTML codes, based on matplitlib.BASE_COLORS and .CSS4_COLORS.
    """
    
    # https://matplotlib.org/examples/color/named_colors.html
    return dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)


def _get_chain_list(traj, chains):
    """
    If chains is "all" returns range for number of chains
    
    Else returns tauren.commons.int_slicer()
    """
    # TypeError: object of type 'list_iterator' has no len()
    number_of_chains = sum(1 for _ in traj.topology.chains)
    
    if chains == "all":
        chain_list = list(range(number_of_chains))
    
    else:
        chain_list = tcommons.int_slicer(chains, number_of_chains)
    
    log.debug(f"<chain_list>: {chain_list}")
    
    return chain_list


def _get_colors(colors):
    """
    Return an interator based on colors.
    """
    if colors is "None":
        colors = it.cycle(mcolors.BASE_COLORS)
    
    else:
        colors = it.cycle(colors.split(","))
        
    return colors
