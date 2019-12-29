"""Library wide core utils."""
import os
from pathlib import Path as _Path


class Path(type(_Path())):
    """
    Extends Python's `Path object`_ interface.

    .. _Path object: https://docs.python.org/3/library/pathlib.html
    """
    def str(self):
        """
        String version of Path.
        
        Returns
        -------
        str
           ``os.fspath(self)``.
        """
        return os.fspath(self)
    
    def myparents(self):
        """
        A list of the path parent folders.
    
        Returns
        -------
        list
            Parent paths. Name file or folder are excluded.
        """
        return self.resolve().parents[0]

