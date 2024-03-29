"""Library wide core utils."""
import os
from functools import wraps
from pathlib import Path as _Path


issues = 'https://taurenmd.readthedocs.io/en/latest/contributing.html#reporting-and-requesting'  # noqa: E501

controlled_exit = 'Your run could not finish.'
CONTACTUS = f'Write us an issue explaining your problem at: {issues}'


class Path(type(_Path())):
    """
    Extends Python's `Path object`_ interface.

    .. _Path object: https://docs.python.org/3/library/pathlib.html
    """

    def str(self):
        """
        Represent path as string.

        Alias to ``os.fspath(self)``.
        
        Returns
        -------
        str
           ``os.fspath(self)``.
        """
        return os.fspath(self)
    
    def myparents(self):
        """
        List of the path parent folders.

        Alias to ``pathlib.Path.resolve().parents[0]``.

        Returns
        -------
        list
            Parent paths. Name file or folder are excluded.
        """
        return self.resolve().parents[0]


def add_reference(ref):
    """
    Add reference decorator.

    Example
    -------

        >>> @add_reference(str)
        >>> def myfunct():
        >>>     ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            references.add(ref)
            result = func(*args, **kwargs)
            return result
        return wrapper
    return decorator


references = set()


ref_taurenmd = "* Cite taurenmd according to: https://taurenmd.readthedocs.io/en/latest/citing.html\n"  # noqa: E501
"""How to cite taurenmd project."""

ref_openmm = "* Data loaded with [OpenMM](http://openmm.org/)"  # noqa: E501
"""Command-line docstring to reference OpenMM package."""

ref_mdt = "* MD data accessed and/or processed using [MDTraj](https://mdtraj.org/)\n"  # noqa: E501
"""Command-line docstring to reference MDTraj package."""

ref_mda = "* MD data accessed using [MDAnalysis](https://www.mdanalysis.org).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis package."""

ref_mda_selection = "* selection commands follow MDAnalysis [selection nomenclature](https://www.mdanalysis.org/docs/documentation_pages/selections.html#).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis selection commands."""

ref_mda_unwrap = "* unwrap performed by MDAnalysis [unwrap](https://www.mdanalysis.org/docs/documentation_pages/core/groups.html?highlight=unwrap#MDAnalysis.core.groups.AtomGroup.unwrap).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis selection.unwrap method."""

ref_mda_alignto = "* align performed by MDAnalysis [unwrap](https://www.mdanalysis.org/docs/documentation_pages/analysis/align.html?highlight=alignto#MDAnalysis.analysis.align.alignto).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis alignto function."""

ref_numpy = '* Matrix operations performed with [Numpy](https://www.scipy.org/citing.html).'  # noqa: E501
"""Command-line docstring to reference numpy lib."""
