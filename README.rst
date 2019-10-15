taurenmd
========

*This is an experimental project.*

Command-line and library interface for analysis routines in Molecular Dynamics.

Installation
============

Currently, installing :code:`taurenmd` does not install its dependencies, you should install them manually, it is better to create a new Python environment for that; :code:`taurenmd` dependencies are described in the :code:`requirements.yml` file.

First clone this repository::

    git clone https://github.com/joaomcteixeira/taurenmd

Inside the repository folder and, if you use `Anaconda`_::

    conda env create -f requirements.yml

Activate the environment::

    conda activate taurenmd

To install :code:`taurenmd` ::
    
    python setup.py develop

To keep your installation up to date::

    git pull

from within the project folder.

Available interfaces
====================

The following interfaces are available, access them through the command-line:

* taurenmd_cif2noHOH
* taurenmd_traj2pdb

Acknowledges
============

The concept of this project is largely inspired in the `pdb-tools`_ *one script one action* idea.
Thanks to `JoaoRodrigues`_ for all the mentoring on MD!

.. _pdb-tools: https://github.com/haddocking/pdb-tools
.. _JoaoRodrigues: https://github.com/JoaoRodrigues
.. _Anaconda: https://www.anaconda.com/distribution/

