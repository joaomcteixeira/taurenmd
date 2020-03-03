taurenmd
========

.. image:: https://raw.githubusercontent.com/joaomcteixeira/taurenmd/master/docs/logo/taurenmd_logo_black_readme.png

.. start-description

.. image:: https://img.shields.io/pypi/v/taurenmd.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/taurenmd

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3551990.svg
    :target: https://doi.org/10.5281/zenodo.3551990
    :alt: Zenodo

.. image:: https://img.shields.io/readthedocs/taurenmd/latest?label=Read%20the%20Docs
    :target: https://taurenmd.readthedocs.io/en/latest/index.html
    :alt: Read the Docs (latest)

.. image:: https://img.shields.io/travis/joaomcteixeira/taurenmd/master?label=TravisCI
    :target: https://travis-ci.org/joaomcteixeira/taurenmd
    :alt: Travis master branch

.. image:: https://codecov.io/gh/joaomcteixeira/taurenmd/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/joaomcteixeira/taurenmd
    :alt: Codecov master branch

.. image:: https://api.codeclimate.com/v1/badges/d69e2e9866338d88955c/test_coverage
   :target: https://codeclimate.com/github/joaomcteixeira/taurenmd/test_coverage
   :alt: Test Coverage

.. image:: https://api.codeclimate.com/v1/badges/d69e2e9866338d88955c/maintainability
   :target: https://codeclimate.com/github/joaomcteixeira/taurenmd
   :alt: Code Climate

.. image:: https://img.shields.io/codeclimate/tech-debt/joaomcteixeira/taurenmd?label=Code%20Climate%20tech%20debt
    :target: https://codeclimate.com/github/joaomcteixeira/taurenmd
    :alt: Code Climate technical debt

.. image:: https://img.shields.io/pypi/wheel/taurenmd.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/taurenmd

.. image:: https://img.shields.io/pypi/pyversions/taurenmd.svg
    :alt: Supported versions
    :target: https://pypi.org/project/taurenmd

.. image:: https://img.shields.io/github/commits-since/joaomcteixeira/taurenmd/v0.8.9/master
    :alt: Commits since latest release
    :target: https://github.com/joaomcteixeira/taurenmd/compare/v0.8.9...master

.. image:: https://img.shields.io/pypi/dm/taurenmd?label=PyPI%20Downloads
    :alt: PyPI - Downloads
    :target: https://pypistats.org/packages/taurenmd

**A command-line interface for analysis routines of Molecular Dynamics data.**
  
**Taurenmd** provides an easy, flexible and extensible, **command-line** interface for the most common *(and not so common)* routines of analysis and representation of Molecular Dynamics (MD) data.

It bridges the gap between the highly complex (and powerful) Python libraries available for analysis of MD data and the *non-developer* users that lack the programming skills to perform a thorough and proficient use those libraries. *But not only*, **taurenmd** also facilitates high throughput operations, even to those proficient *devs*, because complex executions are reduced to single argument-rich command-lines that can be concatenated or aliased.

**Taurenmd** wraps around feature-rich and powerful MD analysis libraries such as `MDAnalysis <https://www.mdanalysis.org/>`_ and `MDTraj <http://mdtraj.org>`_ *(but not only)*, combining them to extract the best of *those worlds*. We use these libraries to access and extract MD data and calculate observables, and we have also added our own routines of analysis when needed. When using this software, you **should** cite **taurenmd** together with the dependencies used, please read our `Citing <https://taurenmd.readthedocs.io/en/latest/citing.html>`_ page for a detailed explanation.

Though designed to perform as a command-line user-directed interface, all **taurenmd**'s core functions are openly distributed and documented. Currently, there are already `several command-line interfaces available <https://taurenmd.readthedocs.io/en/latest/reference/clients.html>`_, some that perform only single tasks, while others allow complex setups, all are *one-liners*.

With this said, **taurenmd** aims to be a flexible and extensible peace of software, built as simple and modular as we can think of, to *agile* the incorporation of new functionalities as needed.

.. end-description

Documentation
=============

**taurenmd** full documentation is available at: https://taurenmd.readthedocs.io, read there:

#. how to install
#. usage
#. citing
#. *etc...*
