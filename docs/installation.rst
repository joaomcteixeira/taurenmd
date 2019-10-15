============
Installation
============

Currently, installing :code:`taurenmd` does not install its dependencies, you should install them manually, it is better to create a new Python environment for that; :code:`taurenmd` dependencies are described in the :code:`requirements.yml` file.

First clone this repository::

    git clone https://github.com/joaomcteixeira/taurenmd

Inside the cloned repository folder and, if you use `Anaconda`_::

    conda env create -f requirements.yml

Activate the environment::

    conda activate taurenmd

To install :code:`taurenmd` ::
    
    python setup.py develop

To keep your installation up to date, just updated the cloned folder::

    git pull


.. _Anaconda: https://www.anaconda.com/distribution/
