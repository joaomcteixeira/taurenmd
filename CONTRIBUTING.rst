============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given. You can contribute
from the scope of an user or as a core Python developer.

Reporting and Requesting
========================

Bug reports
-----------

When `reporting a bug <https://github.com/joaomcteixeira/taurenmd/issues>`_ please use one of the provided issue templates if applicable, otherwise just start a blank issue and describe your situation.

Documentation improvements
--------------------------

**taurenmd** could always use more documentation, whether as part of the
official docs, in docstrings, or even on the web in blog posts,
articles, and such. Write us a `documentation issue <https://github.com/joaomcteixeira/taurenmd/issues/new?assignees=joaomcteixeira&labels=documentation&template=documentation.md&title=%5BDOCUMENTATION%5D>`_ describing what you
would like to see improved in here, and, if you can do
it, just `Pull Request <https://github.com/joaomcteixeira/taurenmd/pulls>`_ your proposed updates ``:-)``.

Feature requests and feedback
-----------------------------

The best way to send feedback is to file an issue using the `feature template <https://github.com/joaomcteixeira/taurenmd/issues/new?assignees=joaomcteixeira&labels=enhancement&template=feature_request.md&title=%5BFEATURE%5D>`_.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions are welcome :)

Code Development
================

General instructions
--------------------

Though *taurenmd* can be used without relying on Anaconda to manage its dependencies,
there are actually some dependencies that are only available in the Anaconda ecosystem.
Therefore, the full development of *taurenmd* is bound to Anaconda until this situation is solved,
you can read further on this in our :ref:`installation page <Installation>`.

Also, by using Conda, our CI integration is much faster because the heavy dependencies
are handled just faster in Anaconda than in pure PyPI.

For these reasons, the development of *taurenmd* is bound to Anaconda. However, it is
possible to make some gymnastics and avoid using Conda and still contribute
to the develop *taurenmd*; if this is your case, please write us an issue and
we will help you out.

From this moment on we will assume you are using Anaconda as your Python package manager.

1. Create a new environment and install *taurenmd* dependencies, **only its dependencies**, following the instructions
for :ref:`Anaconda in our Installation <With Anaconda>` page. **Important**, change the name of the
environment in the ``yml`` file before actually creating the environment, for example,
name it ``taurenmddev``.

    1.1. Remember to activate the new environment before proceeding to the installation.

2. *taurenmd* relies on `tox <https://tox.readthedocs.io/en/latest/>`_ to manage continuous integration (CI) and collaborative development; install it together with `tox-conda <https://github.com/tox-dev/tox-conda>`_::

    conda install -c conda-forge tox
    conda install -c conda-forge tox-conda

3. Fork `taurenmd <https://github.com/joaomcteixeira/taurenmd>`_ (look for the "Fork" button on the top right corner of `our repository <https://github.com/joaomcteixeira/taurenmd>`_).


4. `Clone <https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository>`_ your forked repository to your local machine::

    git clone https://github.com/YOUR-USER-NAME/taurenmd.git <destination folder> 

5. Navigate to the fork folder and create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

6. Install a development version of your development branch, remember to active the ``taurenmddev`` environment beforehand::

    python setup.py develop

Now you can make your changes locally.


7. When you're done making changes run all the checks and docs builder with **tox** one command::

    tox

And correct the errors reported. You can run individual test environment with *tox*, for example,
to test syntax and code style and building related tasks::

    tox -e check 

to test documentation::

    tox -e docs

to perform coverage-reported tests::

    tox -e py37
    # or if you are implementing on python3.6
    tox -e py36

altogether, for example::

    tox -e check,docs

The ``tox`` command  will run testing environments for both ``py36`` and ``py37``,
it will fail if you don't have both interpreters installed. If that is the case,
just ignore the errors for that env.

8. Commit your changes and push your branch to your *taurenmd fork* on GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

9. `Submit a pull request through the GitHub website <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`_.

A New Command-line client
-------------------------

One of the most natural and straightforward ways to contribute to *taurenmd* is
to develop a new command-line interface that implements a new analysis routine;
in this way this routine becomes available in the *taurenmd* catalog.
We provide a command-line client
`template file <https://github.com/joaomcteixeira/taurenmd/blob/master/src/taurenmd/_cli_template.py>`_
from which you can start developing your own new command-line client,
just copy the template file to a new file named ``cli_NAME.py`` and follow
the instructions provided as comments in that same file, the instructions
will guide you throughout all required steps.
Found the template file under ``src/taurenmd/`` folder.

Pull Request Guidelines
-----------------------

If you need some code review or feedback while you're developing the code just make a pull request.

For merging, you should:

1. Make sure your PR passes all ``tox`` tests [1]_.
2. Update documentation when there's new API, functionality etc.
3. Add a note to ``CHANGELOG.rst`` about the changes.
4. Add yourself to ``AUTHORS.rst``.

.. [1] If you don't have all the necessary python versions available locally you can rely on Travis - it will
       `run the tests <https://travis-ci.org/joaomcteixeira/taurenmd/pull_requests>`_ for each change you add in the pull request.

       It will be slower though ...

Continuous Integration
======================

This project follows Continuous Integration (CI) good practices (let us know if something can be improved). As referred in the previous section, CI environment is provided by `tox <https://tox.readthedocs.io/en/latest/>`_ in combination with `tox-conda <https://github.com/tox-dev/tox-conda>`_. All *tox* testing environments run on `Travis-CI <https://travis-ci.org/joaomcteixeira/taurenmd>`_; there, we check for code style integrity, documentation, tests and test coverage. CI configuration is defined in the `tox.ini <https://github.com/joaomcteixeira/taurenmd/blob/master/tox.ini>`_ and in the `.travis.yml <https://github.com/joaomcteixeira/taurenmd/blob/master/.travis.yml>`_ files.

Currently, we do not provide testing for Windows and MacOSX platforms. *taurenmd* depends on several research libraries and we cannot, and should not, attempt to guarantee proper installation of those libraries on all platforms. Therefore we decided to provide full test coverage just for Linux systems where we know those libraries operate fully. You may wish to read our `Installation page for more comments on this matter <https://taurenmd.readthedocs.io/en/latest/installation.html#supported-platforms>`_.
