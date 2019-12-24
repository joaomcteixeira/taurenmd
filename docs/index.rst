Welcome to taurenmd Documentation!
===================================

.. image:: https://img.shields.io/readthedocs/taurenmd/latest?label=Read%20the%20Docs
    :target: https://taurenmd.readthedocs.io/en/latest/index.html
    :alt: Read the Docs (latest)

.. image:: https://img.shields.io/github/commits-since/joaomcteixeira/taurenmd/v0.7.0/develop
    :alt: Commits since latest release
    :target: https://github.com/joaomcteixeira/taurenmd/compare/v0.7.0...develop

**A command-line interface for analysis routines in Molecular Dynamics (MD).**

**Taurenmd** provides an easy, flexible and extensible, interface for the most common (and not so common) routines of analysis and data representation for Molecular Dynamics.

It bridges the gap between highly complex (and powerful) Python libraries for MD data analysis and *non-developer* users that have to struggle to use them for routine needs. Also, :code:`taurenmd` facilitates high throughput operations, even to those proficient *devs*, because complex executions are reduced to a single argument-rich command line which can be easily run under different environments, the same way you would run :code:`ls` in different folders.

**Taurenmd** wraps around feature-rich and powerful MD analysis libraries such as :mda:`/` and :mdtraj:`/` (but not only), combining them to extract the best of *those worlds*. We use these libraries to access and extract MD data and calculate observables, but we have also added our own routines of analysis when needed.

With this said, **taurenmd** aims to be a flexible and extensible peace of software, built as simple and modular as we can think of, to *agile* the incorporation of new functionalities as needed.

**taurenmd** aims to be a command-line user directed interface, but its core functions are openly available and distributed in the library architecture. Currently, there are already several *cmds* available, some that perform only one tasks, others that allow a more complex setup, all *one-liners*.

You can read now through the **contents** bellow.


Contents
========

.. toctree::
   :maxdepth: 1

   installation
   usage
   citing
   reference/index
   contributing
   authors
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


