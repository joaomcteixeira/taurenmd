.. Tauren-MD documentation master file, created by
   sphinx-quickstart on Mon Feb 11 15:49:48 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tauren-MD's documentation!
=====================================

An interface that streamlines analisis routines for Molecular Dynamics.

Designed mostly for non-developer users, but that serves also those
proficient Python programmers.

Tauren-MD wraps around high performance MD analysis libraries,
currently: `MDTraj`_ and `OpenMM`_ (`MDAnalysis`_ to be implemented);
and other Python libraries, such as `matplotlib`_, to allow simple yet
dynamic analysis workflows.

Tauren-MD contains sets of predefined functions for *loading*,
*transforming* and *exporting* trajectories, for calculating restraints
from trajectories and also predifined plotting routines to maximize the
output quality. Tauren-MD is an *À la carte* menu where the user can
easily chose what she/he wants from hers/his trajectories.

Tauren-MD does not require any programming skills to be used, you can
easily configure the ``tauren_config.json`` file with the routines you
wish to execute and run it via ``taurenmd`` executable file (created
after installation).

.. toctree::
   :maxdepth: 2
   :caption: Contents:
    
   docs/installation
   docs/guided_user_interface
   docs/command_line_interface
   docs/modules
   docs/license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _OpenMM: https://github.com/pandegroup/openmm
.. _MDAnalysis: https://www.mdanalysis.org/
.. _matplotlib: https://matplotlib.org/

