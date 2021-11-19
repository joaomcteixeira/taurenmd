Citing
======

.. image:: https://joss.theoj.org/papers/10.21105/joss.02175/status.svg
    :target: https://doi.org/10.21105/joss.02175
    :alt: joss

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3551990.svg
    :target: https://doi.org/10.5281/zenodo.3551990
    :alt: Zenodo

If you use **taurenmd**, please cite it in your publications and research/work environment, referring always to the main publication:

* Teixeira, J. M., (2020). taurenmd: A command-line interface for analysis of Molecular Dynamics simulations. Journal of Open Source Software, 5(50), 2175. https://doi.org/10.21105/joss.02175

**BibText**::

    @article{Teixeira2020,
      doi = {10.21105/joss.02175},
      url = {https://doi.org/10.21105/joss.02175},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {50},
      pages = {2175},
      author = {João M.C. Teixeira},
      title = {taurenmd: A command-line interface for analysis of Molecular Dynamics simulations.},
      journal = {Journal of Open Source Software}
    }

This project is also indexed at `Zenodo <https://doi.org/10.5281/zenodo.3551990>`_ and, if needed, you can complement the above citation in two by:

#. Referring to the whole project::

    João M.C. Teixeira. joaomcteixeira/taurenmd: Zenodo. http://doi.org/10.5281/zenodo.3551990

#. **or** referring the exact version used, for example::

    João M.C. Teixeira. (2019, December 26). joaomcteixeira/taurenmd: v0.7.2 (Version v0.7.2). Zenodo. http://doi.org/10.5281/zenodo.3593004

Citing is one of the best ways to support this project.

Citing Dependencies
-------------------

When using and citing ``taurenmd``, you **SHOULD by all means** cite the Molecular Dynamics (MD) analysis libraries, and others, with which ``taurenmd`` operates to perform the executions you have used. These MD libraries are the :ref:`tauremd software dependencies <How taurenmd manages its dependencies>`. Taurenmd uses different libraries for the different :ref:`client interfaces <Command-line interfaces>`, each command-line interface documentation has a ``References`` section that indicates third party libraries used that you should cite.


After each command execution, the command used for that execution is appended to the ``.taurenmd.cmd`` file (see :ref:`logging <Logging>`). Also a reference to each third pary project used during that execution is appended after the command registry, follow the citing instruction to properly cite the related projects. 

Bellow, links to citing instructions of other research projects **taurenmd** uses as dependencies, alphabetical order:

#. `MDAnalysis citing <https://www.mdanalysis.org/pages/citations/>`_
#. `MDTraj citing <http://mdtraj.org/1.9.3/index.html?highlight=citing#citation-doi-for-citing-mdtraj>`_
#. `matplotlotlib citing <https://matplotlib.org/3.1.1/citing.html>`_
#. `Numpy citing <https://www.scipy.org/citing.html>`_
#. `OpenMM citing <https://simtk.org/projects/openmm>`_
#. `PyQuaternion citing <https://github.com/KieranWynn/pyquaternion>`_
#. `python-bioplottemplates citing <https://github.com/joaomcteixeira/python-bioplottemplates/>`_
