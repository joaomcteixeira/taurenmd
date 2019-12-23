Usage
=====

**Taurenmd** provides a command-line interface to many routines in Molecular Dynamics data analysis. **Taurenmd** uses and wraps around *third-party* libraries to handle MD data, please refer to our :ref:`citing` section.

Our :ref:`Command-line interfaces` page documents in detail each and every of the available `command-line` interfaces; to have access to a short list execute.::

    taurenmd

this command displays a list of the available `taurenmd` interfaces.

Here we make use of the :ref:`trajedit` as an example.::

    taurenmd trajedit -h

With :code:`trajedit` you can edit a trajectory in many different ways. For example, let's convert a trajectory to another format:::

    taurenmd trajedit topology.pdb trajectory.xtc -d new_trajectory.dcd

the above command reads the original `trajectory.xtc` file and outputs the new :code:`new_trajectory.dcd`. You can also use :ref:`trajedit` to reduce the trajectory size, say by slicing every 10 frames:::

    taurenmd trajedit topology.pdb trajetory.xtc -d traj_p10.xtc -p 10 -O

the :code:`-p` option revers to the slicing step size, in this case of :code:`10`. Likewise, you can pass a *start* (:code:`-s`) and an *end* (:code:`-e`) parameters:::

    taurenmd trajedit topology.pdb trajectory.xtc -d traj_s50_e500_p10.xtc -s 50 -e 500 -p 10 -O

in the above examples, the option :code:`-O` disables the creationg a topology file for the new trajectory, because by default a new topology file is created from the first frame of the new trajectory.

You can also use :ref:`trajedit` to extract a specfic frame from a trajectory:::

    taurenmd trajedit topology.pdb trajectory.xtc -d frame40.pdb -s 40 -e 41 -O

but you can instead use the :ref:`fext` interface for this.


