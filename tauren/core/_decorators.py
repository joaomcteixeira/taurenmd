"""
Tauren-MD decorators.

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

import functools
import mdtraj as md
from pathlib import Path


def validate_file_paths(func):
    """
    Validates paths. Paths should exists and be files.
    Raises Erros otherwise.
    
    Used in functions where all positional arguments are paths.
    """
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        
        for file_ in args:
            
            path_ = Path(file_)
            
            if not path_.exists():
                raise FileNotFoundError(f"'{path_}' does NOT exist.")
            
            if not path_.is_file():
                raise ValueError(f"'{path_}' is NOT a file.")
            
        return func(*args, **kwargs)
    
    return wrapper


def validate_trajectory(func):
    """
    Valiates trajectory types.
    """
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        
        if not isinstance(args[0], md.Trajectory):
            raise TypeError("Not a valid trajectory type")
        
        return func(*args, **kwargs)
    
    return wrapper
