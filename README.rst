Tauren-MD
=========

|doi| |version|

.. |doi| image:: https://img.shields.io/badge/DOI-10.5281%2Fzenodo.2579632-informational.svg
   :target: https://zenodo.org/badge/latestdoi/152575798

.. |version| image:: https://img.shields.io/static/v1.svg?label=version&message=0.6.0&color=orange
   :target: https://github.com/joaomcteixeira/Tauren-MD/releases

**Tauren-MD is an interface that streamlines analysis routines for Molecular Dynamics.**

Description
-----------

Tauren-MD was designed to facilitate the usage of scientific MD analysis libraries to non-developer users, though it can be very useful also for the proficient Pythonistas out there.

Tauren-MD wraps around high performance MD analysis libraries, such as: `MDTraj`_, `MDAnalysis`_, `OpenMM`_ (and implementing others...); and it contains its own routines for data representation and export, such as curated plotting templates through `matplotlib`_.

In this way, Tauren-MD attempts to be an *à la carte* menu where the user can easily chose from a predefined list of *actions* what she/he wants from hers/his trajectories: a configuration file can be easily setup with the desired routines and executed via the ``taurenmd`` executable file (created after installation).

.. _version: https://semver.org/#semantic-versioning-200
.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _MDAnalysis: https://www.mdanalysis.org/
.. _OpenMM: https://github.com/pandegroup/openmm
.. _matplotlib: https://matplotlib.org/

Documentation
-------------

The complete and latest Tauren-MD project description and full documentation is available *online* at the `project's website`_.

If you have cloned this repository or downloaded any of the previous `releases`_ you can access the corresponding version documentation *offline*, simply open the ``index.html`` file stored in ``docs`` folder with your favourite web-browser.

.. _`project's website`: https://joaomcteixeira.github.io/Tauren-MD/
.. _releases: https://github.com/joaomcteixeira/Tauren-MD/releases

Citing
------

Cite **Tauren-MD** project as:

- João M.C. Teixeira. joaomcteixeira/Tauren-MD. Zenodo. http://doi.org/10.5281/zenodo.2579632.

Or, if you want to cite a specific version you have used, refer to the `DOI of that version`_, for example:

- João M.C. Teixeira. (2019, February 28). joaomcteixeira/Tauren-MD: v0.5.2 (Version v0.5.2). Zenodo. http://doi.org/10.5281/zenodo.2580076

You **SHOULD** by all means cite the Molecular Dynamics (MD) analysis libraries upon which Tauren-MD wraps. For example, by selecting ``"mdtraj"`` as ``"traj_type"`` in the `configuration file`_ you are making use of it, therefore you should also cite MDTraj along with Tauren-MD; the same applies for any other library Tauren-MD uses and which you have chosen to use.

.. _`DOI of that version`: https://zenodo.org/record/2580076#.XH_8jYVw30o
.. _`configuration file`: https://joaomcteixeira.github.io/Tauren-MD/taurenhtml/html/rstfiles/configuration_file.html#trajectory-type

`Please read further`_ on more details on how to cite Tauren-MD's dependencies.

.. _`Please read further`: https://joaomcteixeira.github.io/Tauren-MD/taurenhtml/html/rstfiles/citing.html

License
-------

The entire Tauren-MD code comes with no liability and is licensed under the `GPL-3.0`_.

.. image:: https://raw.githubusercontent.com/joaomcteixeira/Tauren-MD/master/docs/img/gpl3_logo.png
    :target: https://www.gnu.org/licenses/gpl-3.0.en.html

.. _GPL-3.0: https://github.com/joaomcteixeira/Tauren-MD/blob/master/LICENSE

