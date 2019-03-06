.. _inst:

Installation
============

Automated Installation (Recommended)
------------------------------------

Tauren-MD installation is supported by the :treeoflife:`Tree-of-Life project <>`, which was designed with the very specific :treevision:`vision <>` of facilitating the installation process of scientific software to non-developer users, but it serves perfectly for proficient Pythonistas as well. The installation process is compatible with **Linux**, **MacOS** and **Windows**. To install a stand-alone Taren-MD, run in your ``Terminal`` (command-line):

::

   python install_tauren-md.py

and follow the very simple instructions. If you are on Windows try *double-click*.

This process will install a local :anaconda:`Anaconda <>` distribution and environment with all the Python libraries required to run Tauren-MD. Additionally, two executable files are created in the ``bin`` folder: the ``taurenmd`` and the ``update`` files. If you are on Windows, these files will have the ``.py`` extension.

Read further on :ref:`how to run and update Tauren-MD <run_tauren>`.

Manual installation
-------------------

If you are a proficient Pythonista, either using :pypi"`PyPI <>` or :anaconda:`Anaconda <>`, and want to install Tauren-MD Python library within your Python distribution without using the :treeoflife:`Tree-of-Life project <>`, use the ``install/taurenmd.yml`` file to create the Tauren-MD environment and add the Tauren-MD folder to the ENV `site-packages` or to your system PATH.

However, in this way the ``bin`` executable files won't be created. These scripts are stored in the ``ìnstallation/executables.py`` file as raw strings, copy those to new files to generate your own executable files.

.. todo::
    Add Tauren-MD package to PyPI and/or Conda repository.
