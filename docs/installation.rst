Installation
============

**taurenmd** is written in, and depends on projects written in, `Python <https://www.python.org>`_; therefore, its installation process is based on the Python installation routines and related community-available tools. Find this project's:

#. `package at PyPI <https://pypi.org/project/taurenmd/>`_
#. `GitHub source repository <https://github.com/joaomcteixeira/taurenmd>`_

Dependencies
------------

Installing :code:`taurenmd` does **not** install **all** its dependencies. **Why?** Because taurenmd relies on large and complex libraries required to manage the Molecular Dynamics (MD) data, such as `MDAnalysis <https://www.mdanalysis.org>`_ and `MDTraj <https://mdtraj.org/>`_, and we acknowledge that:

1. it lies outside this project scope to guarantee proper installation of those complex (non-pure Python) libraries, even relying on PyPI or Anaconda.
2. many MD researchers actually do not work on the stable-release version of those *libraries*, instead they work on:

  * cutting edge *development* versions,
  * forked versions,
  * source-compiled versions
  * So installing the deployed version of those libraries by default might be counter productive [1]_
  * For those who want to work with the stable-release version of these dependencies, their respective websites provide very straightforward ways to install them, you should refer to those, *please continue reading*, and just continue with taurenmd installation on top of them.

3. platform compatibility issues (read further)
4. lastly and minor, not all dependencies are required for every *taurenmd command*,

For these reasons, we have decided **not** to install these large dependencies alongside with *taurenmd* installation; other minor dependencies are automatically installed, though. Bellow a list of the dependencies you should install yourself prior to installing *taurenmd*.

.. note::
    
    Please read this whole page before installing the following dependencies.

#. `MDAnalysis Installation instructions <https://www.mdanalysis.org/pages/installation_quick_start/>`_
#. `MDTraj installation instructions <http://mdtraj.org/1.9.3/installation.html>`_
#. `OpenMM installation <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>`_
#. `Numpy <https://numpy.org/>`_, is installed together with the above dependencies, so you should not need to reinstall it again, just stick to the version compatible with the 3 libraries, this should be managed automatically by your Python package manager. Nonetheless, and for your interest, **taurenmd** requires *Numpy* but it is not installed along with the main installation.

Other dependencies installed automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other dependencies that are indeed automatically installed are listed bellow:

#. `python-bioplottemplates <https://github.com/joaomcteixeira/python-bioplottemplates>`_
#. `pyquaterion <http://kieranwynn.github.io/pyquaternion/>`_

Supported Platforms
-------------------

**taurenmd** is designed to run natively under any `platform compatible with Python <https://pythondev.readthedocs.io/platforms.html>`_ (paths are not hard coded ``;-)``). However, the libraries **taurenmd** depends on may or may not be compatible with all OS platforms, and we are not responsible for providing such compatibility or support; you should choose a platform compatible with all the required Molecular Dynamics analysis libraries used by *taurenmd*. We can guarantee **taurenmd** works fully with all its dependencies under Ubuntu 18.04 LTS. You may wish to read further on how continuous integration is :ref:`managed for this project <continuous integration>`.

Installation steps
------------------

From previous environment
~~~~~~~~~~~~~~~~~~~~~~~~~

If you use Molecular Dynamics for your research, odds are you have already the :ref:`above mentioned packages installed <Dependencies>`; if this is the case, you can just install *taurenmd* on top of them, run: ``pip install taurenmd`` on your MD analysis Python environment.

From scratch
~~~~~~~~~~~~

1. Create a new Python environment running Python 3.6 (or 3.7), such procedure varies depending on which Python package manager you use:

  1.1. if your are using `Anaconda`_ as your Python package manager, `read here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_. At the day of writing it would be::
    
    conda create -n taurenmdenv python=3.7

  1.2. if you are using `PyPI`_ as you Python package manager, `read here instead <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_.

3. Activate the newly created environment; if you are with Anaconda::

    conda activate taurenmdenv

4. Install each of the *large library* dependencies, visit their respective website for instructions, see `Dependencies`_ section.

5. Either you use Anaconda or PyPI, install **taurenmd** running the following command::

    pip install taurenmd

6. To update **taurenmd** simply::

    pip install -upgrade taurenmd

7. At this moment you should have all ``taurenmd`` command-line interfaces available on your system, run on your ``terminal``::

    taurenmd

8. In case something is failing, please write us an `Issue <https://github.com/joaomcteixeira/taurenmd/issues>`_ explaining your situation.

From GitHub
```````````

If you are a proficient Pythonista you might want to install **taurenmd** from a development branch on GitHub. If that is the case you might not need to read this section because you  already know well what to do; nonetheless, let's go through it:

.. note::

    ``taurenmd`` follows Semantic Version 2.0, meaning that every single new addition to the master branch gets released on PyPI with a new version number.
    Therefore, installing from the ``master`` GitHub branch actually adds no benefit to installing from PyPI.

#. install the MD analysis libraries as described in the above sections
#. clone our repository: :code:`git clone https://github.com/joaomcteixeira/taurenmd`
#. place yourself in the new :code:`taurenmd` folder, in Linux-like systems: :code:`cd taurenmd`.
#. ``git checkout the-branch-you-want-to-use``
#. install **taurenmd** with the following command: :code:`python setup.py develop`
#. in the future, to keep your installation up to the latest:

  #. pull repository updates from the upstream repository: :code:`git pull` (from within :code:`taurenmd` git folder)
  #. just in case something special was added, repeat :code:`python setup.py develop`

.. [1] though this could be disabled by the user using the ``--no-deps`` option
.. _PyPi: https://pypi.org/
.. _Anaconda: https://www.anaconda.com/distribution/
