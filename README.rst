Tauren-MD
=========

.. image:: https://zenodo.org/badge/152575798.svg
   :target: https://zenodo.org/badge/latestdoi/152575798

`version`_: 0.5.2

**Tauren-MD is an interface that streamlines analisis routines for Molecular Dynamics.**

Description
-----------

Tauren-MD was designed to facilitate the usage of scientific MD analysis libraries to non-developer users, though it can be very useful also for the proficient Pythonistas out there. 

Tauren-MD wraps around high performance MD analysis libraries, such as: `MDTraj`_, `MDAnalysis`_, `OpenMM`_ (and implementing others...); and it contains its own routines for data representation and export, such as curated plotting templates through `matplotlib`_.

In this way, Tauren-MD attempts to be an *à la carte* menu where the user can easily chose from a predifined list of *actions* what she/he wants from hers/his trajectories: a configuration file can be easily setup with the desired routines and executed via the ``taurenmd`` executable file (created after installation).

.. _version: https://semver.org/#semantic-versioning-200
.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _MDAnalysis: https://www.mdanalysis.org/
.. _OpenMM: https://github.com/pandegroup/openmm
.. _matplotlib: https://matplotlib.org/

Documentation
-------------

The complete and latest Tauren-MD project description and full documentation is available `online`_. If you have cloned this repository or downloaded any of the previous `releases_` you can access the corresponding version documentation offline, simply open the ``index.html`` file stored in ``docs`` folder with your favourite web-browser.

.. _online: https://joaomcteixeira.github.io/Tauren-MD/
.. _releases_: https://github.com/joaomcteixeira/Tauren-MD/releases

License
-------

The entire Tauren-MD code comes with no liability and is licensed under the `GPL-3.0`_.

.. image:: https://raw.githubusercontent.com/joaomcteixeira/Tauren-MD/master/docs/img/gpl3_logo.png
    :target: https://www.gnu.org/licenses/gpl-3.0.en.html

.. _GPL-3.0: https://github.com/joaomcteixeira/Tauren-MD/blob/master/LICENSE
