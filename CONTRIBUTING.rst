============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

Bug reports
===========

When `reporting a bug <https://github.com/joaomcteixeira/taurenmd/issues>`_ please include:

    * Your operating system name and version.
    * Any details about your local setup that might be helpful in troubleshooting.
    * Detailed steps to reproduce the bug.

Documentation improvements
==========================

taurenmd could always use more documentation, whether as part of the
official taurenmd docs, in docstrings, or even on the web in blog posts,
articles, and such.

Feature requests and feedback
=============================

The best way to send feedback is to file an issue at https://github.com/joaomcteixeira/taurenmd/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions are welcome :)

Development
===========

To contribute to the development of taurenmd, set up a local environment:

1. Create a new :code:`taurenmd` environment and install all its dependencies but do NOT install *taurenmd* itself, see :ref:`Installation` page.

2. taurenmd relies on `tox <https://tox.readthedocs.io/en/latest/>`_ to manage continuous integration (CI) and collaborative development; install it together with `tox-conda <https://github.com/tox-dev/tox-conda>`_::

    # with Anaconda
    conda install -c conda-forge tox
    conda install -c conda-forge tox-conda
    
    # with PyPi
    pip install tox
    pip install tox-conda

3. Fork `taurenmd <https://github.com/joaomcteixeira/taurenmd>`_ (look for the "Fork" button).

4. Clone your fork locally::

    git clone https://github.com/YOUR-USER-NAME/taurenmd.git <destination folder> 

5. Navigate to the fork folder and create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

6. Install a development version of your development branch::

    python setup.py develop

   Now you can make your changes locally.

7. When you're done making changes run all the checks and docs builder with **tox** one command::

    tox

5. Commit your changes and push your branch to your *taurenmd fork* on GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

6. `Submit a pull request through the GitHub website <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`_.

Pull Request Guidelines
-----------------------

If you need some code review or feedback while you're developing the code just make the pull request.

For merging, you should:

1. Include passing tests (run ``tox``) [1]_.
2. Update documentation when there's new API, functionality etc.
3. Add a note to ``CHANGELOG.rst`` about the changes.
4. Add yourself to ``AUTHORS.rst``.

.. [1] If you don't have all the necessary python versions available locally you can rely on Travis - it will
       `run the tests <https://travis-ci.org/joaomcteixeira/taurenmd/pull_requests>`_ for each change you add in the pull request.

       It will be slower though ...

Tips
----

You can run individual test environment with tox, for example, to test lint::

    tox -e check 

to test documentation::

    tox -e docs

to perform coverage-reported tests::

    tox -e py37-cover
