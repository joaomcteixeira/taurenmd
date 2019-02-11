Run Tauren-MD
=============

Tauren-MD is built such that all its functionalities can be used without
requiring any programming skills from the user.

The only necessary requirement is to configure the
``tauren_config.json`` file according to your preferences and run the
``taurenmd`` file in ``bin`` folder. See bellow how to do this.

taurenmd executable
-------------------

To get help on how to run the ``taurenmd`` file simply run:

::

   python taurenmd -h

   optional arguments:
     -h, --help            show this help message and exit
     -c CONFIG, --config CONFIG
                           Tauren-MD configuration JSON file. If not provided
                           runs the default config: simply loads trajectory and
                           shows information
     -traj TRAJECTORY, --trajectory TRAJECTORY
                           Trajectory file (('.xtc', '.nc', '.trr', '.h5',
                           '.pdb', '.binpos', '.dcd'))
     -top TOPOLOGY, --topology TOPOLOGY
                           Topology file (('.pdb', '.cif'))
     -tt {mdtraj}, --trajtype {mdtraj}
                           Library to use as trajectory type, in lower case.
                           Current libraries available: MDTraj
                           (http://mdtraj.org/1.9.0/)

Tauren-MD configuration file
----------------------------

The Tauren-MD configuration file, `tauren_config.json`_, is a JSON
dictionary. It is composed of two main blocks (dictionaries):

1. "Input Data"
2. "Actions"

Input Data
~~~~~~~~~~

The ``"input_data"`` dictionary inside the ``tauren_config.json`` file,
contains the ``PATHS`` to the necessary input data. In the current
version, those are ``trajectory`` and ``topology`` file paths. These
files should be as described by `tauren.load.load_traj()`_.

If the input data is not provided in the JSON file, it should be
provided via the OPTIONAL ARGUMENTS ``--trajectory`` and ``--topology``
of the ``taurenmd`` executable. If input data is provided in the JSON
file but also in the OPTIONAL ARGUMENTS, the OPTIONAL ARGUMENTS will
prevail.

Actions
~~~~~~~

The ``"actions"`` dictionary contains all the *actions* that will be
performed SEQUENTIALLY by the Tauren-MD workflow ``:-)``.

The ``actions`` dictionary ``keys`` are the name of the action to be
perform, and each key's ``value`` is a dictionary of ``kwargs`` for the
function that is associated with that action. For example:

.. code:: json

   "actions": {
       "remove_solvent": {
           "exclude":null
           },
       
       "reduce_equidistant": {
           "step":10
           },
       
       "frames2PDB": {
           "frames": "all",
           "suffix": "_"
           }
       }

with this configuration, Tauren-MD will sequentially perform:

1. remove solvent molecules from trajectory
2. reduce the trajectory in equidistant frames spaced by 10 frames
3. export all frames from the reduced trajectory

Reorder actions
~~~~~~~~~~~~~~~

Most importantly, you can reorder the actions by simple reordering the
``actions`` dictionary. For example:

.. code:: json

   "actions": {
       "reduce_equidistant": {
           "step":10
           },
       
       "remove_solvent": {
           "exclude":null
           },
       
       "frames2PDB": {
           "frames": "all",
           "suffix": "_"
           }
       }

with this configuration, Tauren-MD will sequentially perform:

1. reduce the trajectory in equidistant frames spaced by 10 frames
2. remove solvent molecules from trajectory
3. export all frames from the reduced trajectory

Add and remove actions
~~~~~~~~~~~~~~~~~~~~~~

To activate (add) or inactivate (remove) an action simply add or remove
that action's dictionary from the main ``actions`` dictionary, for
example:

Repeating actions
~~~~~~~~~~~~~~~~~

It's very easy to repeat an action. Simply, add tailling underscores
``_`` to the action name such as no action name is repeated. For
example, if you want to export the trajectory in different formats:

.. code:: json

   "actions": {
       "reduce_equidistant": {
           "step":10
           },
       
       "remove_solvent": {
           "exclude":null
           },
       
       "frames2PDB": {
           "frames": "all",
           "suffix": "_"
           },
       
       "save_traj": {
           "file_name": "traj_OUTPUT.dcd",
           "overwrite":true
           },
       
       "save_traj_": {
           "file_name": "traj_OUTPUT.dcd",
           "overwrite":true
           }
       }

Other notes
~~~~~~~~~~~

Also, when adding actions, REMEMBER to add the trailing comma ``,``
after the action command, EXCEPT for the last action - see the examples
above.

List of Actions
===============

Bellow the list of all user configurable actions in
Tauren-MD.

All actions are functions inside Tauren-MD modules, that are stored
inside ``tauren`` directory; if you are a **developer** please go
forward exploting it if you wish ``;-)``.

If you are an **user** and just want to use Tauren-MD, the templates
bellow describe the actions available. Just copy the actions dicitionary
to the configuration file file accordingly to your preferences. You may
wish to read first the what was stated above.

Transform
---------

Actions that transform the trajectory.

remove solvent
~~~~~~~~~~~~~~

.. code:: json

   "remove_solvent": {
       "exclude":null
       }

image molecules
~~~~~~~~~~~~~~~

.. code:: json

   "try_image_molecules": {
       "anchor_molecules": null,
       "other_molecules": null,
       "sorted_bonds": null,
       "make_whole": null
       }

reduce equidistant
~~~~~~~~~~~~~~~~~~

.. code:: json

   "reduce_equidistant": {
       "step":10
       },

slice
~~~~~

.. code:: json

   "slice": {
       "start": 1,
       "stop": 100,
       "step": 1
       }

Export
------

Frames to PDB
~~~~~~~~~~~~~

.. code:: json

   "frames2PDB": {
       "frames": "all",
       "prefix": "_"
       }

save trajectory
~~~~~~~~~~~~~~~

.. code:: json

   "save_traj": {
       "file_name": "traj_OUTPUT.dcd",
       "overwrite":true
       }

Data calculation and plotting
-----------------------------

RMSDs of combined chains
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: json

   "produce_rmsds_combined_chains": {
               
       "calc_rmsds_combined_chains": {
           "chains": "all",
           "ref_frame": 0
           },
       
       "export_data": {
           "file_name": null,
           "sep": ","
           },
       
       "plot_rmsd_combined_chains": {
           "label": null,
           "suptitle": "Combined Chain RMSDs",
           "x_label": "Frame Number",
           "y_label": "RMSDs",
           "color": "blue",
           "alpha": 0.7,
           "grid": true,
           "grid_color": "lightgrey",
           "grid_ls": "-",
           "grid_lw": 1,
           "grid_alpha": 0.5,
           "legend": true,
           "legend_fs": 6,
           "legend_loc": 4,
           "fig_name": null
           }
       }

RMSDs of separated chains
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: json

   "produce_rmsds_separated_chains": {
               
       "calc_rmsds_separated_chains" : {
           "chains": "all",
           "ref_frame": 0
           },
       
       "export_data": {
           "file_name": null,
           "sep": ","
           },
       
       "plot_rmsd_chain_per_subplot": {
           "labels": null,
           "suptitle": "RMSDs per chain",
           "x_label": "Frame Number",
           "y_label": "RMSDs",
           "colors": null,
           "alpha": 0.7,
           "grid": true,
           "grid_color": "lightgrey",
           "grid_ls": "-",
           "grid_lw": 1,
           "grid_alpha": 0.5,
           "legend": true,
           "legend_fs": 6,
           "legend_loc": 4,
           "fig_name": null
           },
       
       "plot_rmsd_individual_chains_one_subplot": {
           "labels": null,
           "suptitle": "Chains' RMSDs",
           "x_label": "Frame Number",
           "y_label": "RMSDs",
           "colors": null,
           "alpha": 0.7,
           "grid": true,
           "grid_color": "lightgrey",
           "grid_ls": "-",
           "grid_lw": 1,
           "grid_alpha": 0.5,
           "legend": true,
           "legend_fs": 6,
           "legend_loc": 4,
           "fig_name": null
           }

        }

.. _tauren_config.json: https://github.com/joaomcteixeira/Tauren-MD/blob/master/tauren_config.json
.. _tauren.load.load_traj():
