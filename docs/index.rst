.. Tauren-MD documentation master file, created by
   sphinx-quickstart on Mon Feb 11 15:49:48 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tauren-MD's documentation!
=====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Tauren-MD
=========

*version: 0.4.0*

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

::

   bin/taurenmd -c your_config.json

You can access further functionalities via:

::

   bin/taurenmd -h

For programmers, Tauren-MD libraries can also be imported and used as
such, ``import tauren``.

Installation
------------

Tauren-MD installation is supported by the `Tree-of-Life project`_,
which was designed with a very specific `vision`_ aimed at non-developer
users. To install a stand-alone Taren-MD, simply:

::

   python install_tauren-md.py

Documentation
-------------

-  Read the Tauren-MD documentation `here`_, where all the
   functionalities are explained, see how easy it is to use it! ``:-)``

License
-------

The entire Tauren-MD code comes with no liability and is licensed under
the `GPL-3.0`_.

.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _OpenMM: https://github.com/pandegroup/openmm
.. _MDAnalysis: https://www.mdanalysis.org/
.. _matplotlib: https://matplotlib.org/
.. _Tree-of-Life project: https://github.com/joaomcteixeira/Tree-of-Life
.. _vision: https://github.com/joaomcteixeira/Tree-of-Life/blob/master/VISION.md
.. _here: https://github.com/joaomcteixeira/Tauren-MD/wiki
.. _GPL-3.0: https://github.com/joaomcteixeira/Tauren-MD/blob/master/LICENSE
