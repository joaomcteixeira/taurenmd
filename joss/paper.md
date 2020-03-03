---
title: 'taurenmd: A command-line interface for analysis of Molecular Dynamics simulations.'
tags:
  - Python
  - Molecular Dynamics
  - Structural Biology
  - Proteins
  - DNA
  - RNA
  - Biochemistry
authors:
  - name: João M.C. Teixeira
    orcid: 0000-0002-9113-0622
    affiliation: "1, 2"
affiliations:
  - name: "Previous, Biomolecular NMR Laboratory, Organic Chemistry Section, Inorganic and Organic Chemistry Department, University of Barcelona, Baldiri Reixac 10-12, Barcelona 08028, Spain"
    index: 1
  - name: "Current, Program in Molecular Medicine, Hospital for Sick Children, Toronto, Ontario M5G 0A4, Canada"
    index: 2
date: 04 March 2020
bibliography: paper.bib
---

# Summary

Molecular dynamics (MD) simulations of biological molecules have evolved drastically since its application was first demonstrated four decades ago [@Karplus1977] and, nowadays, simulation of systems comprising millions of atoms is possible due to the latest advances in computation and data storage capacity -- and the scientific community's interest is growing [@Hospital2019]. Academic groups develop most of the MD methods and software for MD data handling and analysis. The MD analysis libraries developed solely for the latter scope nicely address the needs of manipulating raw data and calculating structural parameters [@mda1; @mda2; @mdt; @Romo2014; @Roe2013], each with its advantages and drawbacks inherent to their implementation strategies. This diversity enriches the field with a panoply of strategies that the community can utilize.

The MD analysis software libraries widely distributed and adopted by the community share two main characteristics: 1) they are written in pure Python [@CS-R9526], or provide a Python interface; and 2) they are *libraries*: highly versatile and powerful pieces of software that, however, require advanced scripting and understanding to be operated, even for their basic functionalities. While this is the correct approach to develop flexible computational libraries, it creates a barrier between these software packages and the *non-developer* researcher or high throughput practices, particularly for routine data handling. Therefore, the need has emerged to create a platform that efficiently combines the MD libraries available in the Python universe, taking the most out of each, and implements rapid interfaces for routine usage by both experts and non-experts in the field. In response to that, here is presented **taurenmd** (\autoref{fig:logo}), an easy-to-use and extensible ecosystem of command-line interfaces that facilitates complex operations in Molecular Dynamics data analysis by building on top of powerful Python-based MD analysis libraries.

![taurenmd logo.\label{fig:logo}](../docs/logo/taurenmd_logo_black_readme.png)

# Implementation

**taurenmd** provides highly parametrizable command-line interfaces that automate complex operations of Molecular Dynamics (MD) data handling and analysis in unitary executions that represent conceptual ideas, such as, the manipulation of raw MD data or the calculation of structural parameters (*RMSDs, RMSFs, etc...*). Command-line operations are workflows defined by orchestrated single-operation functions. These single-logic functions are coded in the core of *taurenmd*'s library which facilitates their unit-testing and sharing among all interfaces. Therefore, the design of *taurenmd*'s architecture is simple, yet highly modular, flat, easy to read, and extensible. *taurenmd* serves as a hub of routines in continuous growth, where new operations can be implemented and shared among the community in a defined and documented manner. The *taurenmd* project is hosted at GitHub (https://github.com/joaomcteixeira/taurenmd) and is extensively documented at ReadTheDocs (https://taurenmd.readthedocs.io).

To operate on MD data, *taurenmd* uses third-party MD analysis libraries; currently, it imports MDAnalysis [@mda1; @mda2], MDTraj [@mdt], and OpenMM [@OpenMM], and they are used depending on the requirements of each command-line client. But, the design of the program allows facile incorporation of new dependencies to extend or implement new workflows. Finally, though *taurenmd* focuses on enhanced combination of third party libraries, its design leaves room for the implementation of original analysis routines when needed.

The command-line interface of *taurenmd* is hierarchic, where `taurenmd` is the main entry point and the different interfaces exist as subroutines, for example:

```bash
# help instructions for the main taurenmd entry point
$ taurenmd -h
# an execution example
$ taurenmd [SUBROUTINE] [OPTIONS]
# querying help for a specific subroutine
$ taurenmd report -h
```

At the date of publication, *taurenmd* provides ten different command-line interfaces; all of which, with their arguments and functionalities, are thoroughly described in the project's documentation under the "Command-line interfaces" section. Similarly, all individual functional operations provided are open, fully documented, and can be imported and used by other projects if desired.

To invite community contributions, a client template file is provided with detailed instructions to guide the implementation of new command-line workflows. The building blocks required to build command-line clients are also extensively documented in the `libs/libcli.py` module, they are reusable and new blocks can also be added if needed. New logical operations can be implemented in the library core and used in clients. Complete instructions on how to contribute to the project are provided in the documentation. The project provides extensive Continuous Integration tests and explicit instructions for code style and format to guide developers. *taurenmd* follows Semantic Versioning 2.0 and we favor agile development/deployment instead of periodic releases.

# Installation

**taurenmd** is deployed in the Python ecosystem and is available for direct download at PyPI (https://pypi.org/project/taurenmd/):

```bash
$ pip3 install taurenmd[all]
```

*taurenmd* code uses only Python provided interfaces and is, therefore, compatible with any platform able to execute Python. However, the different Molecular Dynamics analysis libraries imported have very different deployment strategies and this project cannot guarantee those will function in all operating systems; it is, however, guaranteed that *taurenmd* works fully on Ubuntu 18.04 LTS running Anaconda as Python package manager. We advise reading the detailed installation instructions provided in the project's documentation.

# Use cases

The *taurenmd* current version has ten command-line interfaces that execute different analysis or data manipulation routines. Extensive usage examples are provided in the documentation website or by the command:

```bash
$ taurenmd -h
```

Here we show how the `trajedit` interface is used for data manipulation and transformation:

```bash
$ taurenmd trajedit topology.pdb trajectory.xtc -d traj_s50_e500_p10.xtc \
> -s 50 -e 500 -p 10 -l 'segid A'
```

The latter extracts a subtrajectory spanning frames 50 to 500 (exclusive) with a step interval of 10 frames, and only for atoms for the `'segid A'` atom group; in this particular case, we make use of MDAnalysis library [@mda1; @mda2] to handle the data.

# Acknowledgements

The initial concept of this project was largely inspired in the pdb-tools project "one script one action" idea [@pdbtools]. The author deeply thanks João P.G.L.M. Rodrigues (ORCID: 0000-0001-9796-3193) for mentoring regarding MD simulations and data analysis and to Susana Barrera-Vilarmau (ORCID: 0000-0003-4868-6593) for her intensive usage of the program since the very first versions and all the discussions, feedback and suggestions on building a user-friendly interface. The project's repository layout and Continuous Integration setup was based on `cookiecutter-pylibrary` repository [@cc] with final personal modifications by J.M.C.T.

# References
