taurenmd
========

**Command-line and library interface for analysis routines in Molecular Dynamics.**
*This is an experimental project.*

.. start-badges

Stable version
--------------

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3551990.svg
    :target: https://doi.org/10.5281/zenodo.3551990
    :alt: Zenodo

.. image:: https://img.shields.io/travis/joaomcteixeira/taurenmd/master?label=TravisCI
    :target: https://travis-ci.org/joaomcteixeira/taurenmd
    :alt: Travis master branch

.. image:: https://ci.appveyor.com/api/projects/status/v9r2032bry817tjh/branch/master?svg=true 
    :target: https://ci.appveyor.com/project/joaomcteixeira/taurenmd
    :alt: Appveyor master branch

.. image:: https://codecov.io/gh/joaomcteixeira/taurenmd/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/joaomcteixeira/taurenmd
    :alt: Codecov master branch

.. image:: https://img.shields.io/coveralls/github/joaomcteixeira/taurenmd/master?label=COVERALLS&logo=COVERALLS
    :target: https://coveralls.io/github/joaomcteixeira/taurenmd
    :alt: Coveralls master

.. image:: https://img.shields.io/codacy/grade/147029f2635e4e62bf670efdef728c28/master?label=Codacy
    :target: https://app.codacy.com/manual/joaomcteixeira/taurenmd/dashboard
    :alt: Codacy master branch

.. image:: https://img.shields.io/readthedocs/taurenmd/stable?label=Read%20the%20Docs
    :target: https://taurenmd.readthedocs.io/en/stable/index.html
    :alt: Read the Docs (stable)

.. image:: https://img.shields.io/pypi/v/taurenmd.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/taurenmd

.. image:: https://img.shields.io/pypi/wheel/taurenmd.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/taurenmd

.. image:: https://img.shields.io/pypi/pyversions/taurenmd.svg
    :alt: Supported versions
    :target: https://pypi.org/project/taurenmd

.. image:: https://img.shields.io/pypi/implementation/taurenmd.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/taurenmd

The stable version is hosted at the `master branch`_.

Notes
~~~~~

- AppVeyor builds only for py36; see `#1`_
- Travis is not yet configured for OSX; see `#2`_ 

Development Branch
------------------

*where new features are tested and code maybe broken :)*

.. image:: https://img.shields.io/travis/joaomcteixeira/taurenmd/develop?label=TravisCI
    :target: https://travis-ci.org/joaomcteixeira/taurenmd
    :alt: Travis-CI latest branch

.. image:: https://ci.appveyor.com/api/projects/status/v9r2032bry817tjh?svg=true
    :target: https://ci.appveyor.com/project/joaomcteixeira/taurenmd
    :alt: Appveyor-CI latest branch

.. image:: https://codecov.io/gh/joaomcteixeira/taurenmd/branch/develop/graph/badge.svg
    :target: https://codecov.io/gh/joaomcteixeira/taurenmd
    :alt: Codecov latest branch

.. image:: https://img.shields.io/coveralls/github/joaomcteixeira/taurenmd/develop?label=COVERALLS&logo=COVERALLS
    :target: https://coveralls.io/github/joaomcteixeira/taurenmd
    :alt: Coveralls latest

.. image:: https://img.shields.io/codacy/grade/147029f2635e4e62bf670efdef728c28/develop?label=Codacy
    :target: https://app.codacy.com/manual/joaomcteixeira/taurenmd/dashboard
    :alt: Codacy latest grade

.. image:: https://api.codeclimate.com/v1/badges/d69e2e9866338d88955c/maintainability
   :target: https://codeclimate.com/github/joaomcteixeira/taurenmd
   :alt: Code Climate

.. image:: https://img.shields.io/codeclimate/tech-debt/joaomcteixeira/taurenmd?label=Code%20Climate%20tech%20debt
    :target: https://codeclimate.com/github/joaomcteixeira/taurenmd
    :alt: Code Climate technical debt

.. image:: https://img.shields.io/readthedocs/taurenmd/latest?label=Read%20the%20Docs
    :target: https://taurenmd.readthedocs.io/en/latest/index.html
    :alt: Read the Docs (latest)

.. image:: https://img.shields.io/github/commits-since/joaomcteixeira/taurenmd/v0.6.0/develop
    :alt: Commits since latest release
    :target: https://github.com/joaomcteixeira/taurenmd/compare/v0.6.0...develop

The latest development is hosted at the `develop branch`_.

Motivation
==========

Provide an easy interface for the most common (and not so common) routines of analysis and data representation for Molecular Dynamics.

Documentation
=============

Read the documentation: https://taurenmd.readthedocs.io

Dependencies
============

:code:`taurenmd` wraps around high performance Molecular Dynamics analysis libraries, such as: `MDAnalysis`_, `MDTraj`_, `OpenMM`_ (and implementing others...); and it contains its own routines for data representation and export, such as curated plotting templates through `matplotlib`_ and `bioplottemplates`_.


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
.. _bioplottemplates: https://github.com/joaomcteixeira/python-bioplottemplates
.. _cookiecutter-pylibrary: https://github.com/ionelmc/cookiecutter-pylibrary
.. _master branch: https://github.com/joaomcteixeira/taurenmd/tree/master
.. _develop branch: https://github.com/joaomcteixeira/taurenmd/tree/develop
.. _#1: https://github.com/joaomcteixeira/taurenmd/issues/1
.. _#2: https://github.com/joaomcteixeira/taurenmd/issues/2
