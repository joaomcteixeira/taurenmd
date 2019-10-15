taurenmd
========

*This is an experimental project.*

Command-line and library interface for analysis routines in Molecular Dynamics.

Dependencies
============

:code:`taurenmd` wraps around high performance Molecular Dynamics analysis libraries, such as: `MDTraj`_, `MDAnalysis`_, `OpenMM`_ (and implementing others...); and it contains its own routines for data representation and export, such as curated plotting templates through `matplotlib`_.

Installation
============

Currently, installing :code:`taurenmd` does not install its dependencies, you should install them manually, it is better to create a new Python environment for that; :code:`taurenmd` dependencies are described in the :code:`requirements.yml` file.

First clone this repository::

    git clone https://github.com/joaomcteixeira/taurenmd

Inside the cloned repository folder and, if you use `Anaconda`_::

    conda env create -f requirements.yml

Activate the environment::

    conda activate taurenmd

To install :code:`taurenmd` ::
    
    python setup.py develop

To keep your installation up to date, just updated the cloned folder::

    git pull

Available interfaces
====================

The following interfaces are available, access them through the command-line:

* taurenmd_cif2noHOH
* taurenmd_traj2pdb

Citing
======

When using and citing :code:`taurenmd`, you **SHOULD** by all means cite the Molecular Dynamics (MD) analysis libraries upon which :code:`taurenmd` wraps. Please read through each project's documentation to understand how to cite them; these projects are linked in the :ref:`Dependencies` header.

Acknowledges
============

The concept of this project is largely inspired in the `pdb-tools`_ *one script one action* idea.
Thanks to `JoaoRodrigues`_ for all the mentoring on MD!

.. _pdb-tools: https://github.com/haddocking/pdb-tools
.. _JoaoRodrigues: https://github.com/JoaoRodrigues
.. _Anaconda: https://www.anaconda.com/distribution/
.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _MDAnalysis: https://www.mdanalysis.org/
.. _OpenMM: https://github.com/pandegroup/openmm
.. _matplotlib: https://matplotlib.org/
