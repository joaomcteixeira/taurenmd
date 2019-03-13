.. _taurenconfig:

Tauren-MD configuration file
============================

How to use the conf file
------------------------

The **Tauren-MD configuration file** is a :wiki:`JSON <JSON>`
dictionary; in *non-developer* words: is a simple text files that follows a certain syntax (indented blocks of text within ``{}`` curly brackets), you can use the main :taurenfiles:`configuration template file <tauren_config.json>` to examine its syntax and as template for yours.

The configuration file is composed of three main blocks, **input data**, **trajectory type** and **actions**, in the following form:

.. code:: json

    {
        "input_data": {},
        "traj_type": "",
        "actions": {}
        }

The ``"input_data"``, ``"traj_type"`` and ``"actions"`` keyswords are *subdictionaries* that should be populated with the desired configuration values. Read bellow to understand how to populate those subdictionaries.

#. :ref:`Input Data <input_data>`
#. :ref:`Trajectory Type <trajtype>`
#. :ref:`Actions <actions>`

.. _input_data:

Input Data
~~~~~~~~~~

The ``"input_data"`` dictionary inside the ``tauren_config.json`` file,
contains the ``PATHS`` to the necessary input data. In the current
version, those are ``trajectory`` and ``topology`` file paths.

.. code:: json

    "input_data": {
            "trajectory": "",
            "topology": ""
        },

If the input data is not provided in the JSON file, it should be
provided via the OPTIONAL ARGUMENTS ``--trajectory`` and ``--topology``
of the ``taurenmd`` executable. If input data is provided in the JSON
file but also in the OPTIONAL ARGUMENTS, the OPTIONAL ARGUMENTS will
prevail.

.. _trajtype:

Trajectory Type
~~~~~~~~~~~~~~~

Selects which MD analysis library will Tauren-MD handle. Currently available libraries are:

- :mdtraj:`MDTraj <mdtraj>` with the keyword ``"mdtraj"``.
- :mdanalysis:`MDAnalysis <mdanalysis>` with keyword ``"mdanalysis"``.

for example:

.. code:: json
    
    "traj_type": "mdtraj",

.. _actions:

Actions
~~~~~~~

The ``"actions"`` dictionary inside the configuration file contains all the **actions** that will be
performed **sequentially** by the Tauren-MD workflow ``:-)``.

The ``actions`` dictionary ``keys`` are the name of the action to be
performed, and each key's ``value`` is a dictionary of ``kwargs`` for the
function that is associated with that action. For example:

.. code:: json

    "actions": {
        "remove_solvent": {
            "exclude":null
            },
       
        "frame_slice": {
            "start": null,
            "end": null,
            "step":10
            },
        
        "frames2file": {
            "frames": "all",
            "prefix": "_",
            "ext": "pdb"
            }
        }

with this configuration, Tauren-MD will sequentially perform:

1. remove solvent molecules from trajectory
2. reduce the trajectory in equidistant steps spaced by 10 frames
3. export all frames from the reduced trajectory

Reorder actions
^^^^^^^^^^^^^^^

Most importantly, you can reorder the actions by simple reordering the
``actions`` dictionary. For example:

.. code:: json

    "actions": {
        "frame_slice": {
            "start": null,
            "end": null,
            "step":10
            },
       
        "remove_solvent": {
            "exclude":null
            },
       
        "frames2file": {
            "frames": "all",
            "prefix": "_",
            "ext": "pdb"
            }
        }

with this configuration, Tauren-MD will sequentially perform:

1. reduce the trajectory in equidistant steps spaced by 10 frames
2. remove solvent molecules from trajectory
3. export all frames from the reduced trajectory

Add and remove actions
^^^^^^^^^^^^^^^^^^^^^^

To **add** or **remove** *actions* from the Tauren-MD run, simply add or remove
that action's dictionary from the main ``actions`` dictionary, taking the above example, if we don't want to *remove the solvent* anymore, we would have the following configuration:

.. code:: json

    "actions": {
        "frame_slice": {
            "start": null,
            "end": null,
            "step":10
            },
       
        "frames2file": {
            "frames": "all",
            "prefix": "_",
            "ext": "pdb"
            }
        }

Deactivate actions
^^^^^^^^^^^^^^^^^^

.. versionadded:: 0.4.1

Alternatively to removing an action from the configuration file, you can **deactivate** an action by adding the ``#`` before that action's name. This allow quick edition of the config file with easy revert. Taking the above example, to stop removing the solvent from the trajectory, simply:

.. code:: json

    "actions": {
        "frame_slice": {
            "start": null,
            "end": null,
            "step":10
            },
       
        "#remove_solvent": {
            "exclude":null
            },
       
        "frames2file": {
            "frames": "all",
            "prefix": "_",
            "ext": "pdb"
            }
        }

Repeating actions
^^^^^^^^^^^^^^^^^

It's very easy to repeat an action. Simply, add tailling underscores
``_`` to the action name such as no action name is repeated. For
example, if you want to export the trajectory in different formats:

.. code:: json

    "actions": {
        "frame_slice": {
            "start": null,
            "end": null,
            "step":10
            },
       
        "remove_solvent": {
            "exclude":null
            },
       
        "save_traj": {
            "file_name": "traj_OUTPUT.dcd",
            "overwrite":true
            },
       
        "save_traj_": {   # <- note here the trailing underscore "_"
            "file_name": "traj_OUTPUT.dcd",
            "overwrite":true
            }
        }

Other notes
~~~~~~~~~~~

Also, when adding actions, REMEMBER to add the trailing comma ``,``
after the action command, EXCEPT for the last action - see the examples
above.

.. code:: json

    "actions": {
        "frame_slice": {
            "start": null,
            "end": null,
            "step":10
            }, # <- trailing comma
       
        "save_traj": {
            "file_name": "traj_OUTPUT.dcd",
            "overwrite":true
            }
        }


List of Actions
---------------

Bellow the list of all user configurable actions in
Tauren-MD.

All actions are functions inside Tauren-MD modules, that are stored
inside ``tauren`` directory; if you are a **developer** or a **Pythonista**, please go
forward exploting it if you wish ``;-)`` - :ref:`Tauren-MD library Documentation <taurenmodules>`.

If you are an **user** and just want to use Tauren-MD, the templates
bellow describe the actions available. Just copy the dicitionaries bellow
to the configuration file ``actions`` main dictionary accordingly to your preferences (you may
wish to read first the what was stated in the above sections).

How to copy the actions
~~~~~~~~~~~~~~~~~~~~~~~

In the Tauren-MD configuration file copy each action dictionary to the main ``action`` dictionary:

.. code:: json

    "actions": {
        
        # COPYT THE ACTION DICTIONARIES HERE
        # MIND THE DICTIONARY FORMMATING RULES
        # SEE THE default tauren_config.json file
        # for an example
        
        }


Transform
~~~~~~~~~

Actions that transform the trajectory.

remove solvent
^^^^^^^^^^^^^^

Removes solvent from the trajectory.

.. code:: json

   "remove_solvent": {
       "exclude":null
       }

image molecules
^^^^^^^^^^^^^^^

.. code:: json

   "try_image_molecules": {
       "anchor_molecules": null,
       "other_molecules": null,
       "sorted_bonds": null,
       "make_whole": null
       }

frame slice
^^^^^^^^^^^

Slices the trajectory in frames. The subsequent operations will consider only the applied selection. You can unselect the slice by providing ``null`` arguments to the ``frame_slice`` parameters.

The ``frame_slice`` considers always slicing from the original trajectory, so you can simply re-slice from the original trajectory by simply setting up another ``frame_slice`` operation.

* **How to slice:**

from frame 1 to frame 100 (included) in steps of 1.

.. code:: json

   "frame_slice": {
       "start": 1,
       "stop": 100,
       "step": 1
       }

* **reslicing:**

from frame 100 to the end.

.. code:: json

   "frame_slice": {
       "start": 100,
       "stop": null,
       "step": 1
       }

* **undo slicing for subsequent operations:**

.. code:: json

   "frame_slice": {
       "start": null,
       "stop": null,
       "step": null
       }

atom select
^^^^^^^^^^^

``atom_select`` selects a group of atoms from the topology in use. Once an atom group has been selected, the following operations will consider only that selection. For example, if ``chainid 1`` or ``resid 1:100`` are selected, a subsequent :ref:`save_traj` will only save the atoms on that selection.

The ``selector`` code is a string (some text inside double quotes ``""``) that must match the specifications of the MD analysis library you chose to use, you sould refer to that library documentation for details: :mdaselections:`MDAnalysis selections <>`, :mdtselections:`MDTraj selections <>`.

.. code:: json

    "atom_select": {
        "selector": "chainid 1"
        }

align trajectory
^^^^^^^^^^^^^^^^

``align_traj`` is available only when choosing ``mdanalysis`` as ``traj_type``.

.. code:: json

    "align_traj": {
        "weights": "mass",
        "inplace": true,
        "file_name": "aligned_traj.dcd"
        }

Export
~~~~~~

Actions that export information from the trajectory.

Frames to file
^^^^^^^^^^^^^^

Export trajectory frames to files.

.. code:: json

   "frames2file": {
       "frames": "all",
       "prefix": "_",
       "ext": "pdb"
       }

.. _save_traj:

save trajectory
^^^^^^^^^^^^^^^

Saves trajectory to a file.

.. code:: json

   "save_traj": {
       "file_name": "traj_OUTPUT.dcd"
       }

Data calculation and plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RMSDs of combined chains
^^^^^^^^^^^^^^^^^^^^^^^^

Calculates, exports and plots combined RMSDs for a set of chains.

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
^^^^^^^^^^^^^^^^^^^^^^^^^

Calculates, exports and plots RMSDs calculated for single chains.

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

