"""
CONTAINS COMMONS FUNCTIONS

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

import sys

from tauren import tlog

log = tlog.get_log(__name__)


def int_slicer(slicer, top, *args, **kwargs):
    """
    Returns a list of (integer) frames from a slicer string.
    
    Slicer can be of the following formats (INCLUSIVE):
        - "0"
        - "0,1,5,100"
        - "0:"
        - ":500"
        - "10:50"
        - "10:100:2"
    """
    
    if not isinstance(slicer, str):
        raise ValueError(f"<slicer> '{slicer}' should be of string type.")
    
    log.debug(f"<slicer>: {slicer}")
    
    slicer_not_valid = f"<slicer> input not valid: '{slicer}'"
    list_of_slices = 0
    
    possible_inputs = (
        slicer.isdigit(),
        (slicer.find(":") < 0 and slicer.count(",") > 0),
        slicer.endswith(":"),
        slicer.startswith(":"),
        len(slicer.split(":")) == 2,
        slicer.count(":") == 2,
        )
    
    dict_of_actions = {
        0: lambda x: list(int(x)),
        1: lambda x: [int(i) for i in x.split(",")],
        2: lambda x: list(range(int(x[:-1]), top + 1)),
        3: lambda x: list(range(0, int(x[1:]) + 1)),
        4: lambda x: list(range(
            int(x.split(":")[0]),
            int(x.split(":")[-1]) + 1)
            ),
        5: lambda x: list(range(
            int(x.split(":")[0]),
            int(x.split(":")[1]) + 1,
            int(x.split(":")[-1]),
            ))
        }
    
    log.debug(possible_inputs)
    
    try:
        which_action = possible_inputs.index(True)
    
    except ValueError:
        log.exception(slicer_not_valid)
        sys.exit(1)
    
    log.debug(f"<index>: {which_action}")
    
    try:
        list_of_slices = dict_of_actions[which_action](slicer)
    
    except ValueError:
        log.exception(slicer_not_valid)
        sys.exit(1)
    
    except Exception:
        log.exception(slicer_not_valid)
        sys.exit(1)
    
    return list_of_slices
