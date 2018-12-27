# -*- coding: utf-8 -*-
"""
LOGGER MODULE.

THIS FILE WAS ADAPTED FROM TREE-OF-LIFE PROJECT (version 1.1.1 - LGPLv3)
AND MODIFIED ACCORDINGLY TO THE NEEDS OF THE TAUREN-MD PROJECT.

Visit the original Tree-of-Life project at:

https://github.com/joaomcteixeira/Tree-of-Life

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
import logging
import sys


class InstallLogger():

    log_file_name = 'treeoflife.log'
    
    def __init__(self, name, log_file_name=None):
        """
        Manages the logging system. Writes INFO to stdout and log file
        and DEBUG to log file.
        
        Parameters:
        
            - name (str): optimal __name__ var
            
            - log_file_name (str): the name of the log file
        """
        
        if log_file_name:
            InstallLogger.log_file_name = log_file_name
        
        self.log_file = InstallLogger.log_file_name
        self.name = name
        
        return
    
    def gen_logger(self):
        """
        Starts, configures and returns logger.
        """
        logger = logging.getLogger(self.name)
        logger.setLevel(logging.DEBUG)
        
        # create a file handler
        debug_ = logging.FileHandler(self.log_file)
        debug_.setLevel(logging.DEBUG)
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        
        # create a logging format
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - \
%(filename)s:%(name)s:%(funcName)s:%(lineno)d - %(message)s')
        debug_.setFormatter(formatter)
        
        # add the handlers to the logger
        logger.addHandler(debug_)
        logger.addHandler(ch)
        
        return logger
