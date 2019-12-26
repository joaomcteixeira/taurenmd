Usage
=====

**Taurenmd** provides a command-line interface to many routines in Molecular Dynamics data analysis, therefore **taurenmd** runs by executing one-line and argument-rich commands on the :code:`terminal`.

**Taurenmd** uses and wraps around *third-party* libraries to handle MD data, please refer to our :ref:`citing` section.

We have several *command* interfaces already implemented, our :ref:`Command-line interfaces` page documents in detail each and every of them. 

:IMPORTANT: Do not forget to activate the python environment where you installed *taurenmd* in case it is not yet activated. Please read through our :ref:`Installation` page.

To access to the complete list of *taurenmd commands* with a summary information for each, execute::

    taurenmd

Here we make use of the :ref:`trajedit` as an example, lets inspect its options::

    taurenmd trajedit -h

With :code:`trajedit` you can edit a trajectory in many different ways. For example, let's convert a trajectory to another format::

    taurenmd trajedit topology.pdb trajectory.xtc -d new_trajectory.dcd

the above command reads the original `trajectory.xtc` file and outputs the new :code:`new_trajectory.dcd`. You can also use :ref:`trajedit` to reduce the trajectory size, say by slicing every 10 frames::

    taurenmd trajedit topology.pdb trajetory.xtc -d traj_p10.xtc -p 10 -O

the :code:`-p` option refers to the slicing step size, in this case of :code:`10` - reads every 10 frames. Likewise, you can pass a *start* (:code:`-s`) and an *end* (:code:`-e`) parameters::

    taurenmd trajedit topology.pdb trajectory.xtc -d traj_s50_e500_p10.xtc -s 50 -e 500 -p 10 -O

in the above examples, the option :code:`-O` disables the creation of a topology file for the new trajectory, while by default for :code:`trajedit` a new topology file is created from the first frame of the new *edited* trajectory.

You can also use :ref:`trajedit` to extract a specific frame from a trajectory::

    taurenmd trajedit topology.pdb trajectory.xtc -d frame40.pdb -s 40 -e 41 -O

but, for this example, you could instead use the :ref:`fext` interface::

    taurenmd fext topology.pdb trajectory.xtc -f 40 -x .pdb -f frame_

Each an every :code:`taurenmd` sub command is available directly as a main command by appending a :code:`t` to its name, for example::

    taurenmd trajedit
    # equals
    ttrajedit


