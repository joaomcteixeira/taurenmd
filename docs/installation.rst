Installation
============

**taurenmd** is written in, and depends on projects written in, `Python <http://www.python.org>`_; therefore, its installation process is based on the Python installation routines and related community-available tools:

#. `package at PyPi <https://pypi.org/project/taurenmd/>`_
#. `GitHub source repository <https://github.com/joaomcteixeira/taurenmd>`_

Dependencies
------------

Installing :code:`taurenmd` does **not** install **all** its dependencies. **Why?** Because we rely on large and complex dependencies required to manage the Molecular Dynamics (MD) data, such as :mda:`/` and :mdtraj:`/`, and we acknowledge that:

1. a third party software (like *taurenmd*) should not attempt to guarantee proper installation of those complex (non-pure Python) libraries on its own (even via PyPi or Anaconda).
2. many MD users actually do not work on the stable-release version of those *libs*, instead they work on:

  * cutting edge *development* versions,
  * forked versions,
  * source-compiled versions

3. for those who want to work with the stable-release version of these dependencies, their respective websites provide very straightforward ways to install them, you should refer to those.
4. lastly and minor, not all dependencies are required for every *taurenmd command*.

For these reasons, we have decided **not** to install these large dependencies together with *taurenmd*; other minor dependencies are automatically installed, though. Bellow a list of the dependencies you should install yourself prior to installating *taurenmd*.

.. note::
    
    Please read this whole page before installing the above mentioned dependencies.

#. `MDAnalysis Installation instructions <https://www.mdanalysis.org/pages/installation_quick_start/>`_
#. `MDTraj installation instructions <http://mdtraj.org/1.9.3/installation.html>`_
#. `OpenMM installation <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>`_
#. `Numpy <https://numpy.org/>`_, is installed together with the above dependencies, so you should not need to reinstall it again, just stick to the version compatible with the 3 libraries, this should be managed automatically by your Python package manager. Nonetheless, and for your interest, **taurenmd** requires *Numpy* but it is not installed along with the main installation.

Other dependencies installed automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other dependencies automatically installed are listed bellow:

#. `python-bioplottemplates <https://github.com/joaomcteixeira/python-bioplottemplates>`_
#. `pyquaterion <http://kieranwynn.github.io/pyquaternion/>`_

Compatible Platforms
--------------------

**taurenmd** is designed to run natively under any `platform compatible with Python <https://pythondev.readthedocs.io/platforms.html>`_ (paths are not hard coded :code:`;-)`). However, the libraries **taurenmd** depends on may or may not be compatible with all OS platforms, you should look for a platform compatible with all the required Molecular Dynamics analysis libraries. We can guarantee **taurenmd** works fully with all its dependencies under Ubuntu 18.04 LTS.

Installation steps
------------------

From previous environment
~~~~~~~~~~~~~~~~~~~~~~~~~

If you use Molecular Dynamics for your research, odds are you have already the :ref:`above mentioned packages installed <Dependencies>`; if that is the case, you can just install *taurenmd* on top of them: `pip install taurenmd`, in your MD analysis Python environment.

From scratch
~~~~~~~~~~~~

1. Create a new Python environment running Python 3.6 (or 3.7), such procedure varies depending on which Python package manager you use:

  1.1. if your are using `Anaconda`_ as your Python package manager, `read here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_. At the day of writing it would be::
    
    conda create -n taurenmdenv python=3.7

  1.2. if you are using `PyPi`_ as you Python package manager, `read there <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_.

3. Activate the newly created environment, if you are with Anaconda::

    conda activate taurenmdenv

3. Install each of the *large library* dependencies, visit their respective website for instructions, see `Dependencies`_ section.

4. Either you use Anaconda or PyPi, install **taurenmd** running the following command::

    pip install taurenmd

From GitHub
```````````

If you are a proficient Pythonista you might want to install **taurenmd** from its *latest* source on GitHub. If that is the case you might not need to read this section because you know well what to do; nonetheless, let's go through it:

#. install the MD analysis libraries as described in the above sections
#. clone our repository: :code:`git clone https://github.com/joaomcteixeira/taurenmd`
#. place yourself in the new :code:`taurenmd` folder, in Linux-like systems: :code:`cd taurenmd`.
#. install **taurenmd** with the following command: :code:`python setup.py develop`
#. in the future, to keep your installation up to the latests:
  #. pull repository updates from the upstream repository: :code:`git pull` (from within :code:`taurenmd` git folder)
  #. just in case something special was added, repeat :code:`python setup.py develop`


.. _PyPi: https://pypi.org/
.. _Anaconda: https://www.anaconda.com/distribution/
