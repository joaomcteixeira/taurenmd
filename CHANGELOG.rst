Changelog
=========

0.3.1
-----

- topology output now defaults to traj name + `frame0.pdb`
- added .myparents() to Path in __init__

0.3.0 (2019-11-06)
------------------

- Created *develop* branch
- Created client for frame extraction: :code:`cli_fext`
- Added option to disable export of frame0 topology in trajedit

0.2.1 (2019-10-26)
------------------

- dropped py35
- separated lib MDAnalysis from MDTraj
- :code:`libio` concerns only general functions
- improved :code:`imagemol` I/O

0.2.0 (2019-10-26)
------------------

- added :code:`cli_report`

0.1.1 (2019-10-26)
------------------

- corrected libio
- trajectory loads based on MDAnalysis now read and concatenate multiple trajectories.

0.1.0 (2019-10-26)
------------------

- added interfaces:
  - :code:`trajedit`
  - :code:`noSol`
  - :code:`imagemol`
  - :code:`rmsd`
  - :code:`cli template`

0.0.0 (2019-10-15)
------------------

* First release on PyPI.
