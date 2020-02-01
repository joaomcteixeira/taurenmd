Installation
============

**taurenmd** is written in, and depends on projects written in, `Python <https://www.python.org>`_; therefore, its installation process is based on the Python installation routines and related community-available tools. Find *taurenmd*:

#. `package at PyPI <https://pypi.org/project/taurenmd/>`_
#. `GitHub source repository <https://github.com/joaomcteixeira/taurenmd>`_

Supported Platforms
-------------------

**taurenmd** is designed to run natively under any `platform compatible with Python <https://pythondev.readthedocs.io/platforms.html>`_ (paths are not hard coded ``;-)``). However, **the libraries taurenmd depends on may or may not be compatible with all OS platforms**, and we are **not** responsible for providing compatibility or support for such libraries. To be able to exploit all its features, you should choose a platform compatible with all the required Molecular Dynamics analysis libraries used by *taurenmd*. :ref:`At the bottom of this page we have a section that describes taurenmd's dependencies <How taurenmd manages its dependencies>`.

We can **guarantee** *taurenmd* works fully with all its dependencies using Anaconda on Ubuntu 18.04 LTS, and we are *positive* (not sure) it will be the same for any system supporting Anaconda.

Installation steps
------------------

From a previous defined environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you use Molecular Dynamics for your research, odds are you have already installed the :ref:`required dependencies <How taurenmd manages its dependencies>`; if this is the case, you can just install *taurenmd* on top of them, run: ``pip install taurenmd`` in your MD analysis Python environment.

From scratch
~~~~~~~~~~~~

To install *taurenmd* from scratch:

With Anaconda
`````````````

If you use `Anaconda`_ as your Python package manager just do the following on your ``terminal``:

1. Download the *taurenmd* Anaconda environment file from our repository::

    curl https://raw.githubusercontent.com/joaomcteixeira/taurenmd/master/requirements.yml -o taurenmdenv.yml 

If for some reason the above does not work just open the link on your WebBrowser and save the text to a file (or save the file).

2. Create a new Anaconda Python environment to host *taurenmd*::

    conda env create -f taurenmdenv.yml

Where ``taurenmdenv.yml`` is the file downloaded in the previous step.

3. Activate the newly created environment::

    conda activate taurenmd

4. Install *taurenmd* itself::

    pip install taurenmd
    
With PyPI
`````````

If you do not use `Anaconda`_ and you actually rely on `PyPI`_ as your package manager, that is also (almost) perfectly fine.

1. Create a new Python environment if you wish following the `official instructions for your running Python version <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment>`_. We do not provide specific commands for these operations because these change with certain frequency, so it is best to refer to the official sources.

2. Install *taurenmd*::

    python -m pip install --upgrade pip wheel
    pip3 install taurenmd[all]

3. What is the problem with the pure PyPI installation?

*taurenmd* relies on OpenMM to read ``.cif`` topology files, and OpenMM is not deployed on PyPI, you need to `install it through its conda channel <https://anaconda.org/omnia/openmm>`_. Therefore, unless you need to load ``.cif`` files you can use *taurenmd* from a pure PyPI installation. Otherwise, you should follow the :ref:`With Anaconda` instructions. :ref:`May be you want to help us out solving this problem :-) <Contributing>`.

4. You should be good to go

Other Platforms
```````````````

We do not provide support for other distribution platforms such as `HomeBrew <https://brew.sh/>`_ or `Chocolatey <https://chocolatey.org/>`_, but may be you can emulate the steps described above for these systems. Feel welcomed to :ref:`improve this documentation with your insights <Contributing>`!


From GitHub
```````````

If you are a proficient Pythonista you might want to install **taurenmd** from a development branch on GitHub. If that is the case you might not need to read this section because you already know well what to do; nonetheless, let's go through it:

.. note::

    ``taurenmd`` follows :ref:`Semantic Version 2.0 <Versioning>`, meaning that every single new addition to the master branch gets released on PyPI with a new version number. Therefore, installing from the ``master`` GitHub branch actually adds no benefit to installing with ``pip``.

#. Install the MD analysis libraries as described in the above sections
#. Clone our repository: ``git clone https://github.com/joaomcteixeira/taurenmd``
#. Place yourself in the new ``taurenmd`` folder, in Linux-like systems: ``cd taurenmd``.
#. ``git checkout -b the-branch-you-want-to-use``
#. Install **taurenmd** with the following command: ``python setup.py develop``
#. In the future, to keep your installation up to the latest:

  #. pull repository updates from the upstream repository: ``git pull`` (from within ``taurenmd`` git folder)
  #. because taurenmd developments are mostly reflected on new interfaces you need to update those as well: ``python setup.py develop``

Running taurenmd
----------------

After installation you can run *taurenmd* with the following command ``:-)``::

    taurenmd

Please read our :ref:`Usage` page for, *whatelse*, usage instructions and examples.

Upgrade
-------

To upgrade *taurenmd* and all its dependencies to the latest version::

   pip3 install -U --force-reinstall taurenmd

Something failed
----------------

In case something is failing during installation, execution or upgrade, please write us an `Issue <https://github.com/joaomcteixeira/taurenmd/issues>`_ explaining your situation.


How taurenmd manages its dependencies
-------------------------------------

By default, installing ``taurenmd`` does **not** install **all** its dependencies. **Why?** Because *taurenmd* relies on large and complex libraries required to manage the Molecular Dynamics (MD) data, such as `MDAnalysis <https://www.mdanalysis.org>`_ and `MDTraj <https://mdtraj.org/>`_, and installing them automatically might not be the optimal solution for every case, for example:

1. Many MD researchers may actually work on:

  * cutting edge *development* versions,
  * forked versions,
  * source-compiled versions

2. There may be platform compatibility issues (read further),
3. Lastly and minor, not all dependencies are required for every *taurenmd command*,

So installing those libraries by default together with *taurenmd* might be counter productive [1]_.

**Nonetheless**, *taurenmd* does provide an easy way to install this dependencies whenever possible and needed. These details are explained in the :ref:`Installation steps` section above.

The dependencies that are kept separate from the default installation process are listed bellow; here, links point to their respective official installation instructions.

#. `MDAnalysis Installation instructions <https://www.mdanalysis.org/pages/installation_quick_start/>`_
#. `MDTraj installation instructions <http://mdtraj.org/1.9.3/installation.html>`_
#. `OpenMM installation <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>`_
#. `Numpy <https://numpy.org/>`_, is installed together with the above dependencies, so you should not need to reinstall it again, just stick to the version compatible with the 3 libraries, this should be managed automatically by your Python package manager. Nonetheless, and for your interest, **taurenmd** requires *Numpy* but it is not installed along with the main installation.

Other dependencies installed automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other dependencies that are indeed automatically installed alongside with *taurenmd* are listed bellow:

#. `python-bioplottemplates <https://github.com/joaomcteixeira/python-bioplottemplates>`_
#. `pyquaterion <http://kieranwynn.github.io/pyquaternion/>`_

.. [1] Dependency installation could be disabled using the ``--no-deps`` flag of ``pip``, but we decided for the other strategy.
.. _PyPi: https://pypi.org/
.. _Anaconda: https://www.anaconda.com/distribution/
