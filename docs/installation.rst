Installation
============

**taurenmd** is written in, and depends on projects written in, `Python <http://www.python.org>`_, therefore, its installation process is based on the Python installation routines and related community-available tools.

Dependencies
------------

Installing :code:`taurenmd` does *not* install *all* its dependencies. **Why?** **taurenmd** uses large and complex dependencies to manage Molecular Dynamics (MD) data, such as :mda:`/` and :mdtraj:`/`. It is a common practice that installation of such libraries is custumized by each user, depending if the user wants to use cutting edge *development*, specific versions, or source-compiled installations. Therefore, we have decided **not** to install these large dependencies during the installation of *taurenmd*, other minor dependencies are automatically installed. Bellow a list of the dependencies you should install yourself prior to installating *taurenmd*.

#. :mda:`/pages/installation_quick_start/`
#. :mdtraj:`/installation.html`
#. :openmm:`/installation.html`
#. `Numpy <https://numpy.org/>`_, is installed together with the above dependencies, so you should not need to reinstall it again, just stick the the version compatible with the 3 libraries, anaconda takes care of this very well. Nonetheless, for your interest, **taurenmd** required *Numpy* but this is not installed along with the main installation.

Other dependencies installed automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Minor dependencies are automatically installed:

#. `python-bioplottemplates <https://github.com/joaomcteixeira/python-bioplottemplates>`_
#. `pyquaterion <http://kieranwynn.github.io/pyquaternion/>`_

Compatible Platforms
--------------------

**taurenmd** is designed to run natively under any `platform compatible with Python <https://pythondev.readthedocs.io/platforms.html>`_ (no paths are hard coded :code:`;-)`). However, the libraries **taurenmd** depends on may or may not be compatible with all OS platforms, you should look for a platform compatible with all the required Molecular Dynamics analysis libraries; so far our experience, **taurenmd** works best under Linux.

Installation steps
------------------

From previous environment
~~~~~~~~~~~~~~~~~~~~~~~~~

If you use Molecular Dynamics for your research odds are you have already the above mentioned packages installed; if that is the case, you can just install *taurenmd* on top of them. Follow to installation steps.

From scratch
~~~~~~~~~~~~

#. Create a new Python enviroment with Python 3.7. Create an enviroment according to:
  #. if your are using `Anaconda`_ as your package manager, `read here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.
  #. if you are using `PyPi`_ as you package manager, `read here <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_.
#. install each of the *large library* dependencies, visit their respective website for instructions, see above section.
#. either you use Anaconda or PyPi, install **taurenmd** running the following command::

    pip install taurenmd

From GitHub
```````````

If you are a proficient Pythonista you might want to install **taurenmd** from its *latest* source on GitHub. If that is the case you might not need to read this section because you know well what to do; nonetheless, let's go through it:

#. install the MD analysis libraries as described in the above section
#. clone our repository: :code:`git clone https://github.com/joaomcteixeira/taurenmd`
#. place yourself in the new :code:`taurenmd` folder, in Linux-like systems: :code:`cd taurenmd`.
#. install **taurenmd** with the following command: :code:`python setup.py develop`
#. to keep your installation up to the latests:
  #. pull repository updates from the upstream repository: :code:`git pull` (from within taurenmd folder)
  #. just in case something special was added, repeat :code:`python setup.py develop`


.. _PyPi: https://pypi.org/
.. _Anaconda: https://www.anaconda.com/distribution/
