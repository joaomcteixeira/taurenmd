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


ref_taurenmd = "* Cite taurenmd according to: https://taurenmd.readthedocs.io/en/latest/citing.html\n"  # noqa: E501
"""How to cite taurenmd project."""

ref_mdt = "* MD data is accessed and/or processed using [MDTraj](https://mdtraj.org/)\n."  # noqa: E501
"""Command-line docstring to reference MDTraj package."""

ref_mda = "* MD data is accessed using [MDAnalysis](https://www.mdanalysis.org).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis package."""

ref_mda_selection = "* selection commands follow MDAnalysis [selection nomenclature](https://www.mdanalysis.org/docs/documentation_pages/selections.html#).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis selection commands."""

ref_mda_unwrap = "* unwrap performed by MDAnalysis [unwrap](https://www.mdanalysis.org/docs/documentation_pages/core/groups.html?highlight=unwrap#MDAnalysis.core.groups.AtomGroup.unwrap).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis selection.unwrap method."""

ref_mda_alignto = "* align performed by MDAnalysis [unwrap](https://www.mdanalysis.org/docs/documentation_pages/analysis/align.html?highlight=alignto#MDAnalysis.analysis.align.alignto).\n"  # noqa: E501
"""Command-line docstring to reference MDAnalysis alignto function."""

ref_plottemplates_param = "* plotting is performed by [python-bioplottemplates plot param function](https://python-bioplottemplates.readthedocs.io/en/latest/reference/plots.html#bioplottemplates.plots.param.plot).\n"  # noqa: E501
"""Command-line docstring to reference python-bioplottemplates param plot."""

ref_plottemplates_labeldots = "* plotting is performed by [python-bioplottemplates plot labeldots function](https://python-bioplottemplates.readthedocs.io/en/latest/reference/plots.html#bioplottemplates.plots.label_dots.plot).\n"  # noqa: E501
"""Command-line docstring to reference python-bioplottemplates labeldots plot."""  # noqa: E501

ref_pyquaternion = "* Quaterion operations are performed with [pyquaterion](http://kieranwynn.github.io/pyquaternion/).\n"  # noqa: E501
"""Command-line docstring to reference pyquaterion lib."""

ref_numpy = '* Matrix operations were performed with [Numpy](https://www.scipy.org/citing.html).'  # noqa: E501
"""Command-line docstring to reference numpy lib."""
