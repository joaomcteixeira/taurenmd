Changelog
=========

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
