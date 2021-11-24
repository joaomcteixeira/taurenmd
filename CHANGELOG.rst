Versioning
==========

This project follows strictly `Semantic Versioning 2.0 <https://semver.org/#semantic-versioning-200>`_ for version control. 

Upon release of version 1.0.0, all additions to the ``master`` branch are done by PR followed by its respective version increment and release on `PyPI <https://pypi.org/project/taurenmd/>`_.

Upon version 0.8, and before version 1, SV2 major version increments are reflected in the *minor* version number, and minor and patch increments are reflected together in the *patch* version number. Everything else follows SV2 rules, in this way users can track backwards incompatibilities if they happen.

Changelog
=========

v0.10.1 (2021-11-24)
------------------------------------------------------------

* Add support for CIF topologies in libmda using openmm.app.pdbxfile
* removes the install*.yml gitactions

v0.10.0 (2021-11-24)
------------------------------------------------------------

* Updates MDAnalysis to version 2.0.0
* Defines versions for all other dependencies
* Updates CI

v0.9.9 (2021-11-22)
------------------------------------------------------------

* Corrects taking/passing of `insort` argument in some clients.

v0.9.8 (2021-11-22)
------------------------------------------------------------

* Upgraded roll, pitch, and yaw angle calculations with torsion angle strategy
* Removed usage of PyQuaternion

v0.9.7 (2021-11-20)
------------------------------------------------------------

* improved trajectory report (#54)

v0.9.6 (2021-02-22)
------------------------------------------------------------

* Remove README badge for build: not needed.
* Added some comments in the workflow files.
* no changes in code, only CI.

v0.9.5 (2021-02-17)
------------------------------------------------------------

* Upgraded CI to Github Actions according to: https://github.com/joaomcteixeira/python-project-skeleton
* Updated README badges

0.9.4 (2020-06-02)
------------------

* Updates documentation with JOSS citation (PR #49)

0.9.3 (2020-05-25)
------------------

* Improves CONTRIBUTION.rst guidelines (PR #46)

0.9.2 (2020-05-17)
------------------

* Client progression is now represented by a progress bar (PR #44)

0.9.1 (2020-05-15)
------------------

* Improves log in `.taurenmd.cmd` (PR #43)

0.9.0 (2020-05-15)
------------------

* Adds `-i` to every `CLI` interface (PR #42)
* `major` version change because `cli_imagemol` lost backwards compatibility

0.8.14 (2020-05-12)
-------------------

* Updates `tox.ini` file for Continuous Integration (PR #40)

0.8.13 (2020-05-12)
-------------------

* Added support for sequence of trajectories in CLIs that use `MDTraj` (PR #39)

0.8.12 (2020-05-04)
-------------------

* PR: #37
* Installs taurenmd directly with Conda env

0.8.11 (2020-04-03)
-------------------

* PR: #33
* Corrected command representation in .taurenmd.cwd adding quotes when needed

0.8.10 (2020-04-02)
-------------------

* PR: #32
* Corrects incorrect usage of MDAnalysis.analysis.alignto function in ``trajedit``.

0.8.9 (2020-03-03)
------------------

* Changed logos, PR #28

0.8.8 (2020-02-03)
------------------

* PRs: #25 #26 #27
* Added taurenmd logo for readthedocs
* Added tauranmd logo in README
* Added taurenmd repository banner
* Improved details in the documentation
* Removed ``.ci`` folder, unnecessary

0.8.7 (2020-02-02)
------------------

* PR #24
* Added PyPI downloads badge
* Improved installation instructions
* Improved and clarified contributing instructions

0.8.6 (2020-01-20)
------------------

* Restructured pip deps: install_requires only takes bioplottemplates and pyquaternion
* two extras_require: `sup` and `md` and `all` which consider both

0.8.5 (2020-01-20)
------------------

* PR #22
* organized dependencies for PyPI
* PyPI only dependencies are referred as install_requires
* MDAnalysis and MDTraj referred in extras_require
* OpenMM left out from pip, only available in Anaconda

0.8.4 (2020-01-19)
------------------

* PR #15
* Added simtk lib import check for controlled failure 
* added error message output for user

0.8.3 (2020-01-19)
------------------

* PR #16 and #19
* corrected argparse autodoc in ReadTheDocs (mock strategy)
* improved tox configuration with better env separation
* #19 reports a communication error between TravisCI and coverage servers

0.8.2 (2020-01-17)
------------------

* Improved CI workflow
  * Dropped COVERALLS
  * Dropped Codacy
  * Setup test-coverage in CodeClimate
  * created `.codeclimate.yml` with explicit configuration
* updated badges

0.8.1 (2020-01-15)
------------------

* PR #14
* Corrected version display in documentation

0.8.0 (2020-01-15)
------------------

* PR #13
* Code architecture improvements
* Complete project main documentation
* Complete library documentation
* command line documented
* Code clean

0.7.2 (2019-12-25)
------------------

* bridged from 0.7.1
* Dropped Appveyor and EXPLICIT Windows support because of #1.
* restructured project GitHub layout. Deprecated develop branch.
* Readthedocs documentation improvements in structure and content.

0.7.0 (2019-12-23)
------------------

* implemented :code:`cli_rotations`, calculates roll, pitch and yaw
    rotation angles of selection.

0.6.0 (2019-12-15)
------------------

* implemented :code:`cli_rmsf` to calculate RMSFs.

0.5.1 (skipped to 0.6.0)
------------------------

* added sort numbered trajs to :code:`cli_trajedit`
* added sort numbered trajectory paths in lib
* improved :code:`cli_imagemol` readability
* added selection in :code:`cli_noSol`

0.5.0 (2019-11-24)
------------------

* created :code:`cli_angle`. Calculates angles between a plane along the trajectory. Plane is given by the three centre_of_geometries of three selections.
* args to plot passed as list are transformed to tuple
* added distance calc and plot interface :code:`cli_distances`
* :code:`trajedit` now saves topology unwrapped

0.4.1 (2019-11-21)
------------------

* renumbered version to 0.4.1. from 0.3.1
* RMSD Cli now calculates for several selections
* Parse plot vars now registers floats
* corrected fext cli entry point
* added align option to trajedit
* topology model writen from first frame of time slicing
* added unwrap() molecule method from MDAnalysis in :code:`trajedit` with respective options
* topology output now defaults to traj name + :code:`frame0.pdb`
* added .myparents() to Path in :code:`__init__`

0.3.0 (2019-11-06)
------------------

* Created *develop* branch
* Created client for frame extraction: :code:`cli_fext`
* Added option to disable export of frame0 topology in trajedit

0.2.1 (2019-10-26)
------------------

* dropped py35
* separated lib MDAnalysis from MDTraj
* :code:`libio` concerns only general functions
* improved :code:`imagemol` I/O

0.2.0 (2019-10-26)
------------------

* added :code:`cli_report`

0.1.1 (2019-10-26)
------------------

* corrected libio
* trajectory loads based on MDAnalysis now read and concatenate multiple trajectories.

0.1.0 (2019-10-26)
------------------

* added interfaces:
  * :code:`trajedit`
  * :code:`noSol`
  * :code:`imagemol`
  * :code:`rmsd`
  * :code:`cli template`

0.0.0 (2019-10-15)
------------------

* First release on PyPI.
