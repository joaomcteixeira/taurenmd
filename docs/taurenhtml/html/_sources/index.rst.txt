.. Tauren-MD documentation master file, created by
   sphinx-quickstart on Mon Feb 11 15:49:48 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*This is an experimental project under strong development, expect small to drastic changes and updates in its structure, architecture and interfaces.*

Welcome to Tauren-MD's documentation!
=====================================

Tauren-MD is an interface that streamlines analisis routines for Molecular Dynamics (MD).

Tauren-MD was designed to facilitate the usage of scientific MD analysis libraries to non-developer users, though it can be very useful also for the proficient Pythonistas out there. 

Tauren-MD wraps around high performance MD analysis libraries, such as: :mdtraj:`MDTraj <mdtraj>`, :pandegroup:`OpenMM <openmm>` (and implementing :prody:`ProDy <>` and others); and it contains its own routines for data representation and export, such as curated plotting templates through :matplotlib:`matplotlib <>`.

In this way, Tauren-MD attempts to be an *à la carte* menu where the user can easily chose from a predifined list of _actions_ what she/he wants from hers/his trajectories: a configuration file can be easily setup with the desired routines and executed via the `taurenmd` executable file (created after installation).

.. toctree::
   :maxdepth: 1
   :caption: Contents:
    
   rstfiles/vision
   rstfiles/installation
   rstfiles/run_taurenmd
   rstfiles/configuration_file
   rstfiles/modules
   rstfiles/license

Special Contents
~~~~~~~~~~~~~~~~
* :ref:`complete_index`
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* :taurengit:`Github repository <>`
