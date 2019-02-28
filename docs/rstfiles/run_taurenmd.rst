.. _run_tauren:
Run Tauren-MD
=============

Tauren-MD is built such that all its functionalities can be used without
requiring any programming skills from the user.

There are two executable files that can be used to run Tauren-MD: ``taurenmd`` and ``taurengui``; on Windows, these files will have the ``.py`` extension. As the name suggests, the former runs command-line while the latter runs the Tauren-MD Guided User Interface.

The only necessary requirement is to configure the
``tauren_config.json`` file according to your preferences and run the
``taurenmd`` file in ``bin`` folder, see bellow how to do this.

Command-line interface
----------------------

To run Tauren-MD command line you should use the ``bin/taurenmd`` script created after installation.

In Linux or MacOS
~~~~~~~~~~~~~~~~~

.. code:: bash

    $ bin/taurenmd
    

On Windows
~~~~~~~~~~

.. code:: bash

    $ bin/taurenmd.py

**NOTE:** do NOT run with ``python bin/taurenmd.py``, because you will be using the main system's Python interpreter which does not have the Tauren-MD dependencies installed.

Examples
~~~~~~~~

Use the ``bin/taurenmd -h`` option to get updated help on the script usage:

:-c CONFIG, --config CONFIG:
    Tauren-MD configuration JSON file. If not provided
    runs the default config: simply loads trajectory and
    shows information

:-traj TRAJECTORY, --trajectory TRAJECTORY:
    Trajectory file (('.xtc', '.nc', '.trr', '.h5',
    '.pdb', '.binpos', '.dcd'))
    
:-top TOPOLOGY, --topology TOPOLOGY:
    Topology file (('.pdb', '.cif'))

:-tt {mdtraj}, --trajtype {mdtraj}:
    Library to use as trajectory type, in lower case.
    Current libraries available: MDTraj
    (http://mdtraj.org/1.9.0/)

``taurenmd`` can be used to run a Tauren-MD configuration file, you can use it to repeat a previous execution or to run a new one from a configured conf file. Read :ref:`here <taurenconfig>` how to configure a Tauren-MD configuration file.

``-traj`` and ``-top`` arguments can be given in-line or a path to those files can be added to the configuration file as well.

.. code:: bash

    $ bin/taurenmd -c <PATH TO CONFIG FILE> -traj <PATH TO TRAJECTORY FILE> -top <PATH TO TOPOLOGY FILE>

Guided User Interface (GUI)
---------------------------

The GUI is under development and is not yet available. Sorry :-(

.. _tauren_config.json: https://github.com/joaomcteixeira/Tauren-MD/blob/master/tauren_config.json
.. _tauren.load.load_traj():
