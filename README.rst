taurenmd
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |appveyor|
        | |coveralls| |codecov|
        | |codacy| |codeclimate|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/taurenmd/badge/?style=flat
    :target: https://readthedocs.org/projects/taurenmd
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/joaomcteixeira/taurenmd.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/joaomcteixeira/taurenmd

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/github/joaomcteixeira/taurenmd?branch=master&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/joaomcteixeira/taurenmd

.. |coveralls| image:: https://coveralls.io/repos/joaomcteixeira/taurenmd/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/joaomcteixeira/taurenmd

.. |codecov| image:: https://codecov.io/github/joaomcteixeira/taurenmd/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/joaomcteixeira/taurenmd

.. |codacy| image:: https://img.shields.io/codacy/grade/147029f2635e4e62bf670efdef728c28.svg
    :target: https://www.codacy.com/app/joaomcteixeira/taurenmd
    :alt: Codacy Code Quality Status

.. |codeclimate| image:: https://codeclimate.com/github/joaomcteixeira/taurenmd/badges/gpa.svg
   :target: https://codeclimate.com/github/joaomcteixeira/taurenmd
   :alt: CodeClimate Quality Status

.. |version| image:: https://img.shields.io/pypi/v/taurenmd.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/taurenmd

.. |wheel| image:: https://img.shields.io/pypi/wheel/taurenmd.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/taurenmd

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/taurenmd.svg
    :alt: Supported versions
    :target: https://pypi.org/project/taurenmd

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/taurenmd.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/taurenmd

.. |commits-since| image:: https://img.shields.io/github/commits-since/joaomcteixeira/taurenmd/v0.1.0.svg
    :alt: Commits since latest release
    :target: https://github.com/joaomcteixeira/taurenmd/compare/v0.1.0...master


.. end-badges

**Command-line and library interface for analysis routines in Molecular Dynamics.**
*This is an experimental project.*

Documentation
=============

Read the documentation: https://taurenmd.readthedocs.io

Dependencies
============

:code:`taurenmd` wraps around high performance Molecular Dynamics analysis libraries, such as: `MDTraj`_, `MDAnalysis`_, `OpenMM`_ (and implementing others...); and it contains its own routines for data representation and export, such as curated plotting templates through `matplotlib`_.


Citing
======

When using and citing :code:`taurenmd`, you **SHOULD** by all means cite the Molecular Dynamics (MD) analysis libraries upon which :code:`taurenmd` wraps. Please read through each project's documentation to understand how to cite them; these projects are linked in the `Dependencies` header.

Acknowledges
============

The concept of this project is largely inspired in the `pdb-tools`_ *one script one action* idea.
Thanks to `JoaoRodrigues`_ for all the mentoring on MD! CI in this repository provided by `cookiecutter-pylibrary`_ with final setup by me.

.. _pdb-tools: https://github.com/haddocking/pdb-tools
.. _JoaoRodrigues: https://github.com/JoaoRodrigues
.. _MDTraj: https://github.com/mdtraj/mdtraj
.. _MDAnalysis: https://www.mdanalysis.org/
.. _OpenMM: https://github.com/openmm/openmm
.. _matplotlib: https://matplotlib.org/
.. _cookiecutter-pylibrary: https://github.com/ionelmc/cookiecutter-pylibrary
