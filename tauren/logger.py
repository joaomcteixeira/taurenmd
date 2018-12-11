"""
MODULE PROVIDES LOGGER TO TAUREN-MD

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
import logging

log_file_name = 'tauren-md.log'


def get_log(name):
    """
    Configures logger for Tauren MD.
    """
    
    log = logging.getLogger(name)
    log.setLevel(logging.DEBUG)
    
    # create a file handler
    debug_ = logging.FileHandler(log_file_name)
    debug_.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    
    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - \
%(filename)s:%(name)s:%(funcName)s:%(lineno)d - %(message)s')
    debug_.setFormatter(formatter)
    
    # add the handlers to the logger
    log.addHandler(debug_)
    log.addHandler(ch)
    
    return log
