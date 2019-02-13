.. Tauren-MD documentation master file, created by
   sphinx-quickstart on Mon Feb 11 15:49:48 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*This is an experimental project under strong development, expect small to drastic changes and updates in its structure, architecture and interfaces.*

Welcome to Tauren-MD's documentation!
=====================================

Tauren-MD is an interface that streamlines analisis routines for Molecular Dynamics.

It is designed mostly for non-developer users, but can serve also those
proficient Python programmers.

Tauren-MD contains sets of predefined functions for *loading*,
*transforming* and *exporting* trajectories from different sources, for *calculating* and *extracting* restraints from trajectories and also predifined plotting routines to maximize the
output quality. Tauren-MD is an *À la carte* menu from which the user can
easily chose what she/he wants from hers/his trajectories by applying the different predifined actions.

Tauren-MD does not require any programming skills to be used, you can
easily configure the ``tauren_config.json`` file with the routines you
wish to execute and run it via ``taurenmd`` executable file (created
after installation).

Tauren-MD does not implement new core functionalities on MD analisis by itself, instead, it wraps around high performance MD analysis libraries that have been developed by the scientific community worldwide, thus providing a simple yet dynamic preconfigured interface to use them. The current version incorporates:

- `MDTraj`_
- `OpenMM`_

but we are looking forward to include also other MD analysis libraries as new functionalities are implemented, such as:

- `MDAnalysis`_
- `ProDy`_
- `nglview`_
- and others...

... other Python libraries are used, such as `matplotlib`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
    
   rstfiles/vision
   rstfiles/installation
   rstfiles/guided_user_interface
   rstfiles/command_line_interface
   rstfiles/modules
   rstfiles/license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _OpenMM: https://github.com/pandegroup/openmm
.. _MDAnalysis: https://www.mdanalysis.org/
.. _Prody: http://prody.csb.pitt.edu/index.html
.. _nglview: https://github.com/arose/nglview
.. _matplotlib: https://matplotlib.org/

